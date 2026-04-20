#include "det4mc.h"
#include <iostream>
#include "ConsoleColorizer.h"
#include "TH2.h"
#include "TFile.h"
#include "plot4mc.h"
#include <optional>
#include "TLorentzVector.h"

static const int eoev = -1;//end of event value (eoev), printed between events in .det file (-1 will separate entries in mass text file)

const std::pair<int,int> det4mc::offsets[] = {
		{112,40},	//detector0 {ringOffset,wedgeOffset}
		{96,32},	//detector1 {ringOffset,wedgeOffset}
		{80,16},	//detector2 {ringOffset,wedgeOffset}
		{64,24},	//detector3 {ringOffset,wedgeOffset}
		{48,0}		//detector4 {ringOffset,wedgeOffset}
	};

det4mc::det4mc(SABRE_Array* SABRE_Array, SPS_Aperture* SPS_Aperture,
			   std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels,
			   TargetEnergyLoss* targetLoss_par1, TargetEnergyLoss* targetLoss_par2, TargetEnergyLoss* targetLoss_par3, TargetEnergyLoss* targetLoss_par4,
			   SABRE_DeadLayerModel* deadLayerLoss_par1, SABRE_DeadLayerModel* deadLayerLoss_par2, SABRE_DeadLayerModel* deadLayerLoss_par3, SABRE_DeadLayerModel* deadLayerLoss_par4,
			   Beamspot* beamspot,
			   TargetAngularStraggler* straggler_par1, TargetAngularStraggler* straggler_par2, TargetAngularStraggler* straggler_par3, TargetAngularStraggler* straggler_par4)
	: SABRE_Array_(SABRE_Array), SPS_Aperture_(SPS_Aperture),
	  SABREARRAY_EnergyResolutionModels_(SABREARRAY_EnergyResolutionModels),
	  targetLoss_par1_(targetLoss_par1), targetLoss_par2_(targetLoss_par2), targetLoss_par3_(targetLoss_par3), targetLoss_par4_(targetLoss_par4),
	  deadLayerLoss_par1_(deadLayerLoss_par1), deadLayerLoss_par2_(deadLayerLoss_par2), deadLayerLoss_par3_(deadLayerLoss_par3), deadLayerLoss_par4_(deadLayerLoss_par4),
	  straggler_par1_(straggler_par1),straggler_par2_(straggler_par2),straggler_par3_(straggler_par3),straggler_par4_(straggler_par4),
	  nevents_(0),
	  hitej_(0), hit1_(0), hit2_(0), hit3_(0),
	  hitBoth23_(0), hitOnlyEj_(0), hitOnly1_(0), hitOnly2_(0), hitOnly3_(0),
	  hitOnly12_(0), hitOnly23_(0), hitOnly13_(0), hitOnly123_(0),
	  onePartHits_(0), twoPartHits_(0), threePartHits_(0), fourPartHits_(0),
	  detectorHits_(SABRE_Array->size()),
	  beamspot_(beamspot)
	{

	}

bool det4mc::ProcessParticle(Particle& p, const Vec3& origin, EventRecorder* rec, plot4mc* plotter, SimConfig* config, std::ostringstream& ss){

	double dth = p.strag->Sample();
	double dph = p.strag->SamplePhi();

	Vec3 traj;
	traj.SetVectorSpherical(1, p.theta*DEGRAD, p.phi*DEGRAD);

	if(config->GetStraggleEnabled(p.id)){
		Vec3 etheta, ephi;
		etheta.SetVectorCartesian(std::cos(p.theta*DEGRAD)*std::cos(p.phi*DEGRAD), std::cos(p.theta*DEGRAD)*std::sin(p.phi*DEGRAD), -std::sin(p.theta*DEGRAD));
		ephi.SetVectorCartesian(-std::sin(p.phi*DEGRAD), std::cos(p.phi*DEGRAD), 0.);

		Vec3 adjustedTrajectory;
		adjustedTrajectory = std::cos(dth*DEGRAD)*traj + std::sin(dth*DEGRAD)*(std::cos(dph*DEGRAD)*etheta + std::sin(dph*DEGRAD)*ephi);
		double adjth = adjustedTrajectory.GetTheta()*RADDEG;
		double adjph = adjustedTrajectory.GetPhi()*RADDEG;
		if(adjph < 0.) adjph += 360.;
		plotter->FillStraggleHistos(p.theta, p.phi, adjth, adjph, dth, dph);
		traj = adjustedTrajectory.Unit();
	} else {
		plotter->FillStraggleHistos(p.theta, p.phi, p.theta, p.phi, 0, 0);
	}

	double theta = traj.GetTheta()*RADDEG;
	double phi = traj.GetPhi()*RADDEG;

	for(size_t i=0; i<SABRE_Array_->size(); i++){

		//check if traj intersects a detector (returns {ring, wedge})
		std::pair<int, int> hit_rw = SABRE_Array_->at(i)->GetOffsetTrajectoryRingWedge(theta*DEGRAD, phi*DEGRAD, origin);

		if(hit_rw.first != -1 && hit_rw.second != -1){

			//calculate energy loss in target
			double e_afterTarget = p.tLoss->ApplyEnergyLoss(p.energy, theta);

			//calculate energy loss in SABRE dead layer
			Vec3 normal = SABRE_Array_->at(i)->GetNormTilted().Unit();
			double e_afterDeadLayer = p.dLoss->ApplyEnergyLoss(e_afterTarget, traj, normal);

			//apply energy resolution and check threshold
			if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit_rw.first, e_afterDeadLayer, p.smearedERing) &&
				(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit_rw.second, e_afterDeadLayer, p.smearedEWedge))){

				auto angles = SABRE_Array_->GetDetectorThetaPhi(offsets[i].first + hit_rw.first, offsets[i].second + hit_rw.second);

				if(angles.has_value()){
					p.theta_meas = angles->first;
					p.phi_meas = angles->second;
					p.detected = true;

					Vec3 localCoords = SABRE_Array_->at(i)->GetHitCoordinatesRandomWiggle(hit_rw.first, hit_rw.second);

					//add to event recorder
					Hit h = {
							 p.id, (int)i, offsets[i].first+hit_rw.first, offsets[i].second+hit_rw.second,
							 hit_rw.first, hit_rw.second, p.smearedERing, p.smearedEWedge, p.theta_meas,
							 p.phi_meas, localCoords.GetX(), localCoords.GetY()
							};
					rec->AddHit(h);

					//write to string stream for .det file output:
					ss << (p.id*100) + i << "\t" << hit_rw.first << "\t" << hit_rw.second << "\t" << p.smearedERing << "\t" <<
					p.smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << "\n";

					detectorHits_[i]++;
					return true;
				}
			}
		}
	}

	//if no detector in the array was triggered:
	p.detected = false;
	p.theta_meas = -666.;
	p.phi_meas = -666.;
	return false;
}

void det4mc::Run(std::ifstream& infile, std::ofstream& outfile, EventRecorder* EventRecorder, plot4mc* RootPlotter, SimConfig* config){

	//initialize particle array
	Particle particles[4] = {
		{1,0,0,0,false,0,0,0,0,targetLoss_par1_,deadLayerLoss_par1_,straggler_par1_},
		{2,0,0,0,false,0,0,0,0,targetLoss_par2_,deadLayerLoss_par2_,straggler_par2_},
		{3,0,0,0,false,0,0,0,0,targetLoss_par3_,deadLayerLoss_par3_,straggler_par3_},
		{4,0,0,0,false,0,0,0,0,targetLoss_par4_,deadLayerLoss_par4_,straggler_par4_}
	};

	while(infile >> particles[0].energy >> particles[0].theta >> particles[0].phi
				 >> particles[1].energy >> particles[1].theta >> particles[1].phi
				 >> particles[2].energy >> particles[2].theta >> particles[2].phi
				 >> particles[3].energy >> particles[3].theta >> particles[3].phi){

		nevents_ += 1;
		if(nevents_ % 50000 == 0) ConsoleColorizer::PrintBlue("Processed " + std::to_string(nevents_) + " events...\n");

		Vec3 reactionOrigin = beamspot_->GeneratePoint();
		EventRecorder->SetReactionOrigin(reactionOrigin.GetX(),reactionOrigin.GetY(),reactionOrigin.GetZ());

		std::ostringstream ss;
		uint8_t hitMask = 0; //empty mask
		int totalHits = 0;

		//begin by explicitly checking particles[0] in SPS aperture
		bool EjInSPS = false;
		double SPS_E, SPS_Theta, SPS_Phi, SPS_RecoilEx;
		Vec3 ejTraj;
		ejTraj.SetVectorSpherical(1.0, particles[0].theta*DEGRAD, particles[0].phi*DEGRAD);
		if(SPS_Aperture_->IsDetected(ejTraj, reactionOrigin)){

			EjInSPS = true;

			SPS_E = SPS_Aperture_->GetSmearedEnergy(particles[0].energy);
			SPS_Theta = SPS_Aperture_->GetSmearedTheta(particles[0].theta);
			SPS_Phi = SPS_Aperture_->GetSmearedPhi(particles[0].phi);

			TLorentzVector beam(0., 0., std::sqrt(2*config->GetBeam().massMeV*config->GetBeamEnergy()), config->GetBeam().massMeV*config->GetBeamEnergy());
			TLorentzVector target(0., 0., 0., config->GetTarget().massMeV);
			double Pej = std::sqrt(2*config->GetEjectile().massMeV*SPS_E);
			TLorentzVector ejectile(Pej*std::sin(SPS_Theta*DEGRAD)*std::cos(SPS_Phi*DEGRAD),
									Pej*std::sin(SPS_Theta*DEGRAD)*std::sin(SPS_Phi*DEGRAD),
									Pej*std::cos(SPS_Theta*DEGRAD),
									config->GetEjectile().massMeV + SPS_E);

			TLorentzVector recoil = beam + target - ejectile;
			SPS_RecoilEx = recoil.M() - config->GetRecoil().massMeV;

		} else {

			EjInSPS = false;

			SPS_E = -666.;
			SPS_Theta = -666.;
			SPS_Phi = -666.;

			SPS_RecoilEx = -666.;

		}

		EventRecorder->UpdateSPS(EjInSPS,SPS_E,SPS_Theta,SPS_Phi,SPS_RecoilEx);


		for(int i=0; i<4; i++){
			particles[i].detected = false;
			EventRecorder->SetKinematics(i, particles[i].energy, particles[i].theta, particles[i].phi);

			ss << particles[i].energy << "\t" << particles[i].theta << "\t" << particles[i].phi;
			if(i != 3) ss << "\t";
			else ss << "\n";

			if(ProcessParticle(particles[i], reactionOrigin, EventRecorder, RootPlotter, config, ss)){
				hitMask |= (1 << i);//set bit corresponding to particle index
				totalHits += 1;
			}
		}

		if(hitMask & 1) hitej_++;
		if(hitMask & 2) hit1_++;
		if(hitMask & 4) hit2_++;
		if(hitMask & 8) hit3_++;
		if((hitMask & 12) == 12) hitBoth23_++;

		switch(totalHits){
			case 1: onePartHits_++; break;
			case 2: twoPartHits_++; break;
			case 3: threePartHits_++; break;
			case 4: fourPartHits_++; break;
		}

		if(hitMask == 1) hitOnlyEj_++;		//0001
		if(hitMask == 2) hitOnly1_++;		//0010
		if(hitMask == 4) hitOnly2_++;		//0100
		if(hitMask == 8) hitOnly3_++;		//1000

		if(hitMask == 6) hitOnly12_++;		//0110
		if(hitMask == 10) hitOnly13_++;		//1010
		if(hitMask == 12) hitOnly23_++;		//1100

		if(hitMask == 14) hitOnly123_++;	//1110

		//finalize event
		ss << eoev;
		outfile << ss.str() << "\n";
		RootPlotter->ProcessTXTOutput(ss.str());
		EventRecorder->FillEvent();
	}
}

long det4mc::GetNumEvents() const {
	return nevents_;
}

long det4mc::GetHitEj() const {
	return hitej_;
}

long det4mc::GetHit1() const {
	return hit1_;
}

long det4mc::GetHit2() const {
	return hit2_;
}

long det4mc::GetHit3() const {
	return hit3_;
}

long det4mc::GetHitBoth23() const {
	return hitBoth23_;
}

long det4mc::GetHitOnlyEj() const {
	return hitOnlyEj_;
}

long det4mc::GetHitOnly1() const {
	return hitOnly1_;
}

long det4mc::GetHitOnly2() const {
	return hitOnly2_;
}

long det4mc::GetHitOnly3() const {
	return hitOnly3_;
}

long det4mc::GetHitOnly12() const {
	return hitOnly12_;
}

long det4mc::GetHitOnly23() const {
	return hitOnly23_;
}

long det4mc::GetHitOnly13() const {
	return hitOnly13_;
}

long det4mc::GetHitOnly123() const {
	return hitOnly123_;
}

long det4mc::GetOnePartHits() const {
	return onePartHits_;
}

long det4mc::GetTwoPartHits() const {
	return twoPartHits_;
}

long det4mc::GetThreePartHits() const {
	return threePartHits_;
}

long det4mc::GetFourPartHits() const {
	return fourPartHits_;
}

const std::vector<long>& det4mc::GetDetectorHits() const {
	return detectorHits_;
}
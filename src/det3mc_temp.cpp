#include "det3mc_temp.h"
#include <iostream>
#include "ConsoleColorizer.h"
#include "TLorentzVector.h"

static const int eoev = -1;//end of event value (eoev), printed between events in .det file (-1 will separate entries in mass text file)

const std::pair<int, int> det3mc_temp::offsets[] = {
	{112,40}, {96,32}, {80,16}, {64,24}, {48,0}
};

	det3mc_temp::det3mc_temp(SABRE_Array* SABRE_Array_, SPS_Aperture* SPS_Aperture_, 
				std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels,
				TargetEnergyLoss* targetLoss_par1,
				TargetEnergyLoss* targetLoss_par2,
				TargetEnergyLoss* targetLoss_par3,
				TargetEnergyLoss* targetLoss_par4,
				SABRE_DeadLayerModel* deadLayerLoss_par1,
				SABRE_DeadLayerModel* deadLayerLoss_par2,
				SABRE_DeadLayerModel* deadLayerLoss_par3,
				SABRE_DeadLayerModel* deadLayerLoss_par4,
				Beamspot* beamspot,
				TargetAngularStraggler* straggler_par1,
				TargetAngularStraggler* straggler_par2,
				TargetAngularStraggler* straggler_par3,
				TargetAngularStraggler* straggler_par4)
	: SABRE_Array_(SABRE_Array_),
	  SPS_Aperture_(SPS_Aperture_),
	  SABREARRAY_EnergyResolutionModels_(SABREARRAY_EnergyResolutionModels),
	  targetLoss_par1_(targetLoss_par1),
	  targetLoss_par2_(targetLoss_par2),
	  targetLoss_par3_(targetLoss_par3),
	  targetLoss_par4_(targetLoss_par4),
	  deadLayerLoss_par1_(deadLayerLoss_par1),
	  deadLayerLoss_par2_(deadLayerLoss_par2),
	  deadLayerLoss_par3_(deadLayerLoss_par3),
	  deadLayerLoss_par4_(deadLayerLoss_par4),
	  straggler_par1_(straggler_par1),
	  straggler_par2_(straggler_par2),
	  straggler_par3_(straggler_par3),
	  straggler_par4_(straggler_par4),
	  nkinevents_(0),
	  nevents_(0),
	  hitejSPS_(0),
	  nohitejSPS_(0),
	  hit1_(0), hit2_(0), hit3_(0),
	  hitBoth23_(0), hitOnly1_(0), hitOnly2_(0), hitOnly3_(0),
	  hitOnly12_(0), hitOnly23_(0), hitOnly13_(0), hitOnly123_(0),
	  onePartHits_(0), twoPartHits_(0), threePartHits_(0),
	  detectorHits_(SABRE_Array_->size()),
	  beamspot_(beamspot)
	{

	}

bool det3mc_temp::ProcessParticle(Particle& p, const Vec3& origin, EventRecorder* rec, plot3mc* plotter, SimConfig* config, std::ostringstream& ss){

	Vec3 traj;
	traj.SetVectorSpherical(1, p.theta*DEGRAD, p.phi*DEGRAD);

	if(config->GetStraggleEnabled(p.id)){
		double dth = p.strag->Sample();
		double dph = p.strag->SamplePhi();

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
	if(phi<0.) phi += 360.;

	for(size_t i=0; i<SABRE_Array_->size(); i++){

		std::pair<int, int> hit_rw = SABRE_Array_->at(i)->GetOffsetTrajectoryRingWedge(theta*DEGRAD, phi*DEGRAD, origin);

		if(hit_rw.first != -1 && hit_rw.second != -1){

			double e_afterTarget = p.tLoss->ApplyEnergyLoss(p.energy, theta);

			Vec3 normal = SABRE_Array_->at(i)->GetNormTilted().Unit();
			double e_afterDeadLayer = p.dLoss->ApplyEnergyLoss(e_afterTarget, traj, normal);

			if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit_rw.first, e_afterDeadLayer, p.smearedERing) && 
				(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit_rw.second, e_afterDeadLayer, p.smearedEWedge))){

				auto angles = SABRE_Array_->GetDetectorThetaPhi(offsets[i].first + hit_rw.first, offsets[i].second + hit_rw.second);

				if(angles.has_value()){
					p.theta_meas = angles->first;
					p.phi_meas = angles->second;
					p.detected = true;

					Vec3 localCoords = SABRE_Array_->at(i)->GetHitCoordinatesRandomWiggle(hit_rw.first, hit_rw.second);

					Hit h = {p.id, (int)i, offsets[i].first+hit_rw.first, offsets[i].second+hit_rw.second,
							 hit_rw.first, hit_rw.second, p.smearedERing, p.smearedEWedge,
							 p.theta_meas, p.phi_meas, localCoords.GetX(), localCoords.GetY()};
							 rec->AddHit(h);

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

void det3mc_temp::Run(std::ifstream& infile, std::ofstream& outfile, EventRecorder* EventRecorder, plot3mc* RootPlotter, SimConfig* config){

	//particles[0] = ejectile, particles[1] = recoil (NOT DIRECTLY MEASURED), particles[2,3] = resonance decay particles
	Particle particles[4] = {
		{1, 0, 0, 0, false, 0, 0, 0, 0, targetLoss_par1_, deadLayerLoss_par1_, straggler_par1_},
		{2, 0, 0, 0, false, 0, 0, 0, 0, targetLoss_par2_, deadLayerLoss_par2_, straggler_par2_},
		{3, 0, 0, 0, false, 0, 0, 0, 0, targetLoss_par3_, deadLayerLoss_par3_, straggler_par3_},
		{4, 0, 0, 0, false, 0, 0, 0, 0, targetLoss_par4_, deadLayerLoss_par4_, straggler_par4_}
	};

	while(infile >> particles[0].energy >> particles[0].theta >> particles[0].phi >> particles[1].energy >> particles[1].theta >> particles[1].phi >> particles[2].energy >> particles[2].theta >> particles[2].phi >> particles[3].energy >> particles[3].theta >> particles[3].phi){

		nkinevents_ += 1;
		if(nkinevents_ % 50000 == 0) ConsoleColorizer::PrintBlue("Processed " + std::to_string(nkinevents_) + " events...\n");

		Vec3 reactionOrigin = beamspot_->GeneratePoint();
		EventRecorder->SetReactionOrigin(reactionOrigin.GetX(), reactionOrigin.GetY(), reactionOrigin.GetZ());

		std::ostringstream ss;
		uint8_t hitMask = 0;
		int totalHits = 0;

		bool EjInSPS = false;
		double SPS_E, SPS_Theta, SPS_Phi, SPS_RecoilEx;
		Vec3 ejTraj;
		ejTraj.SetVectorSpherical(1.0, particles[0].theta*DEGRAD, particles[0].phi*DEGRAD);
		if(SPS_Aperture_->IsDetected(ejTraj, reactionOrigin)){
			EjInSPS = true;

			SPS_E = SPS_Aperture_->GetSmearedEnergy(particles[0].energy);
			SPS_Theta = SPS_Aperture_->GetSmearedTheta(particles[0].theta);
			SPS_Phi = SPS_Aperture_->GetSmearedPhi(particles[0].phi);

			TLorentzVector beam(0., 0., std::sqrt(2*config->GetBeam().massMeV*config->GetBeamEnergy()), config->GetBeam().massMeV+config->GetBeamEnergy());
			TLorentzVector target(0., 0., 0., config->GetTarget().massMeV);
			double Pej = std::sqrt(2*config->GetEjectile().massMeV*SPS_E);
			TLorentzVector ejectile(Pej*std::sin(SPS_Theta*DEGRAD)*std::cos(SPS_Phi*DEGRAD),
									Pej*std::sin(SPS_Theta*DEGRAD)*std::sin(SPS_Phi*DEGRAD),
									Pej*std::cos(SPS_Theta*DEGRAD),
									config->GetEjectile().massMeV + SPS_E);

			TLorentzVector recoil = beam + target - ejectile;
			SPS_RecoilEx = recoil.M() - config->GetRecoil().massMeV;

			hitejSPS_ += 1;

		} else {

			EjInSPS = false;

			SPS_E = -666.;
			SPS_Theta = -666.;
			SPS_Phi = -666.;

			SPS_RecoilEx = -666.;

		}

		EventRecorder->UpdateSPS(EjInSPS, SPS_E, SPS_Theta, SPS_Phi, SPS_RecoilEx);

		for(int i=0; i<4; i++){
			particles[i].detected = false;
			EventRecorder->SetKinematics(i, particles[i].energy, particles[i].theta, particles[i].phi);

			ss << particles[i].energy << "\t" << particles[i].theta << "\t" << particles[i].phi;
			if(i!=3) ss << "\t";
			else ss << "\n";
		}

		if(config->GetSPSCoincidence()){

			if(EjInSPS){

				for(int i=2; i<4; i++){
					if(ProcessParticle(particles[i], reactionOrigin, EventRecorder, RootPlotter, config, ss)){
						int analysisShift = i - 1;
						hitMask |= (1 << analysisShift);
						totalHits += 1;
					}
				}

			} else {
				nohitejSPS_++;
			}

		} else {

			for(int i=0; i<4; i++){
				if(i==1) continue;
				if(i==0 && EjInSPS) continue;

				if(ProcessParticle(particles[i], reactionOrigin, EventRecorder, RootPlotter, config ,ss)){
					int analysisShift = (i == 0) ? 0 : (i -1);
					hitMask |= (1 << analysisShift);
					totalHits += 1;
				}
			}

		}

		if(totalHits > 0) nevents_ += 1;

		if(hitMask & 1) hit1_ += 1;
		if(hitMask & 2) hit2_ += 1;
		if(hitMask & 4) hit3_ += 1;

		if(totalHits == 1) onePartHits_ += 1;
		if(totalHits == 2) twoPartHits_ += 1;
		if(totalHits == 3) threePartHits_ += 1;

		// switch(hitMask){
		// case 4: hitOnly2_ += 1; break;
		// case 8: hitOnly3_ += 1; break;
		// case 12: hitBoth23_ += 1; break;
		// default: break;
		// }
		switch(hitMask){
		case 1: hitOnly1_ += 1; break;
		case 2: hitOnly2_ += 1; break;
		case 4: hitOnly3_ += 1; break;
		case 3: hitOnly12_ += 1; break;
		case 5: hitOnly13_ += 1; break;
		case 6: hitOnly23_ += 1; hitBoth23_ += 1; break;
		case 7: hitOnly123_ += 1; hitBoth23_ += 1; break;
		default: break;
		}

		double sumSABREEnergy = 0.;
		if(hitMask & 1) sumSABREEnergy += particles[0].smearedERing;
		if(hitMask & 2) sumSABREEnergy += particles[2].smearedERing;
		if(hitMask & 4) sumSABREEnergy += particles[3].smearedERing;

		if(totalHits >= 1){
			RootPlotter->FillSumEnergyHisto(hitMask, totalHits, SPS_RecoilEx, sumSABREEnergy);
		}

		ss << eoev;
		// outfile << ss.str() << "\n";
		// RootPlotter->ProcessTXTOutput(ss.str());
		// EventRecorder->FillEvent();

		if(config->GetSPSCoincidence() && EjInSPS){
			outfile << ss.str() << "\n";
			RootPlotter->ProcessTXTOutput(ss.str());
			EventRecorder->FillEvent();
		} else if(config->GetSPSCoincidence() && !EjInSPS){
			EventRecorder->ResetEvent();
		} else if(!config->GetSPSCoincidence()){
			outfile << ss.str() << "\n";
			RootPlotter->ProcessTXTOutput(ss.str());
			EventRecorder->FillEvent();
		}

	}

}
#include "det4mc.h"
#include <iostream>
#include "ConsoleColorizer.h"
#include "TH2.h"
#include "TFile.h"
#include "plot4mc.h"

static const int eoev = -1;//end of event value (eoev), printed between events in .det file (-1 will separate entries in mass text file)

const std::pair<int,int> det4mc::offsets[] = {
		{112,40},	//detector0 {ringOffset,wedgeOffset}
		{96,32},	//detector1 {ringOffset,wedgeOffset}
		{80,16},	//detector2 {ringOffset,wedgeOffset}
		{64,24},	//detector3 {ringOffset,wedgeOffset}
		{48,0}		//detector4 {ringOffset,wedgeOffset}
	};

det4mc::det4mc(std::vector<SABRE_Detector*>& SABRE_Array,
			   std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels,
			   TargetEnergyLoss* targetLoss_par1, TargetEnergyLoss* targetLoss_par2, TargetEnergyLoss* targetLoss_par3, TargetEnergyLoss* targetLoss_par4,
			   SABRE_DeadLayerModel* deadLayerLoss_par1, SABRE_DeadLayerModel* deadLayerLoss_par2, SABRE_DeadLayerModel* deadLayerLoss_par3, SABRE_DeadLayerModel* deadLayerLoss_par4,
			   Beamspot* beamspot,
			   TargetAngularStraggler* straggler_par1, TargetAngularStraggler* straggler_par2, TargetAngularStraggler* straggler_par3, TargetAngularStraggler* straggler_par4)
	: SABRE_Array_(SABRE_Array),
	  SABREARRAY_EnergyResolutionModels_(SABREARRAY_EnergyResolutionModels),
	  targetLoss_par1_(targetLoss_par1), targetLoss_par2_(targetLoss_par2), targetLoss_par3_(targetLoss_par3), targetLoss_par4_(targetLoss_par4),
	  deadLayerLoss_par1_(deadLayerLoss_par1), deadLayerLoss_par2_(deadLayerLoss_par2), deadLayerLoss_par3_(deadLayerLoss_par3), deadLayerLoss_par4_(deadLayerLoss_par4),
	  straggler_par1_(straggler_par1),straggler_par2_(straggler_par2),straggler_par3_(straggler_par3),straggler_par4_(straggler_par4),
	  nevents_(0),
	  hitej_(0), hit1_(0), hit2_(0), hit3_(0),
	  hitBoth23_(0), hitOnlyEj_(0), hitOnly1_(0), hitOnly2_(0), hitOnly3_(0),
	  hitOnly12_(0), hitOnly23_(0), hitOnly13_(0), hitOnly123_(0),
	  onePartHits_(0), twoPartHits_(0), threePartHits_(0), fourPartHits_(0),
	  detectorHits_(SABRE_Array.size()),
	  beamspot_(beamspot)
	{

	}

void det4mc::Run(std::ifstream& infile, std::ofstream& outfile, RootWriter* RootWriter, plot4mc* RootPlotter, bool targetStraggle1, bool targetStraggle2, bool targetStraggle3, bool targetStraggle4){

	// std::cerr << "SABRE_Array_ size: " << SABRE_Array_.size()
	// 		  << ", SABREARRAY_EnergyResolutionModels_ size: " << SABREARRAY_EnergyResolutionModels_.size()
	// 		  << std::endl;


	double e1, theta1, phi1, e2, theta2, phi2, e3, theta3, phi3, e4, theta4, phi4;

	//TH2D *hBeamSpot = new TH2D("hBeamSpot","BeamSpot",200, -0.05, 0.05, 200, -0.05, 0.05);

	while(infile >> e1 >> theta1 >> phi1 >> e2 >> theta2 >> phi2 >> e3 >> theta3 >> phi3 >> e4 >> theta4 >> phi4){

		nevents_ += 1;
		bool detectedEj = false;	//kin4mc particle 1 == ejectile
		bool detected1 = false;		//kin4mc particle 2 == breakup 1 (first decay particle)
		bool detected2 = false;		//kin4mc particle 3 == breakup 2 (second decay particle)
		bool detected3 = false;		//kin4mc particle 4 == breakup 3 (third decay particle -- left over)

		//get reaction origin based on beamspot:
		Vec3 reactionOrigin = beamspot_->GeneratePoint();
		//hBeamSpot->Fill(reactionOrigin.GetX(),reactionOrigin.GetY());
		//Vec3 reactionOrigin = {0.,0.,0.};

		std::ostringstream ss;
		ss << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\t" << e3 << "\t" << theta3 << "\t" << phi3 << "\t" << e4 << "\t" << theta4 << "\t" << phi4 << std::endl;

		if(nevents_%50000==0) ConsoleColorizer::PrintBlue("Processed " + std::to_string(nevents_) + " events...\n");

		RootWriter->SetKinematics(0,e1,theta1,phi1);
		RootWriter->SetKinematics(1,e2,theta2,phi2);
		RootWriter->SetKinematics(2,e3,theta3,phi3);
		RootWriter->SetKinematics(3,e4,theta4,phi4);

		RootWriter->SetReactionOrigin(reactionOrigin.GetX(), reactionOrigin.GetY(), reactionOrigin.GetZ());


		for(size_t i=0; i<SABRE_Array_.size(); i++){
			/*///////////////////////////////////////////////
			//  check if particle 1 (ejectile) is detected //
			///////////////////////////////////////////////*/
			double smearedERing = 0., smearedEWedge = 0.;

			double dtheta1 = straggler_par1_->Sample();
			double dphi1 = straggler_par1_->SamplePhi();

			//original kinematic trajectory
			Vec3 originalTrajectory1;
			originalTrajectory1.SetVectorSpherical(1, theta1*DEGRAD, phi1*DEGRAD);

			//define new basis vectors
			Vec3 etheta1, ephi1;
			etheta1.SetVectorCartesian(std::cos(theta1*DEGRAD)*std::cos(phi1*DEGRAD), std::cos(theta1*DEGRAD)*std::sin(phi1*DEGRAD), -std::sin(theta1*DEGRAD));
			ephi1.SetVectorCartesian(-std::sin(phi1*DEGRAD), std::cos(phi1*DEGRAD), 0.);

			//adjusted trajectory:
			Vec3 adjustedTrajectory1;
			adjustedTrajectory1 = std::cos(dtheta1*DEGRAD)*originalTrajectory1 + std::sin(dtheta1*DEGRAD)*(std::cos(dphi1*DEGRAD)*etheta1 + std::sin(dphi1*DEGRAD)*ephi1);
			adjustedTrajectory1 = adjustedTrajectory1.Unit();

			//now convert back
			double theta1_prime, phi1_prime;
			if(targetStraggle1){

				theta1_prime = adjustedTrajectory1.GetTheta()*RADDEG;
				phi1_prime = adjustedTrajectory1.GetPhi()*RADDEG;
				if(phi1_prime < 0) phi1_prime += 360.;
				RootPlotter->FillStraggleHistos(theta1, phi1, theta1_prime, phi1_prime, dtheta1, dphi1);

			} else {

				theta1_prime = theta1;
				phi1_prime = phi1;
				RootPlotter->FillStraggleHistos(theta1, phi1, theta1_prime, phi1_prime, 0, 0);

			}

			std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta1_prime*DEGRAD,phi1_prime*DEGRAD,reactionOrigin);

			if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detectedEj){
				//std::cout << "here1" << std::endl;
				//apply target energy loss to e1:
				double e1_aftertarget = targetLoss_par1_->ApplyEnergyLoss(e1, theta1_prime);

				//apply dead layer energy loss to e1_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta1_prime*DEGRAD,phi1_prime*DEGRAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e1_afterDeadLayer = deadLayerLoss_par1_->ApplyEnergyLoss(e1_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit1_rw.first, e1_afterDeadLayer, smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit1_rw.second, e1_afterDeadLayer, smearedEWedge)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
					ss << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detectedEj = true;
					hitej_ += 1;
					detectorHits_[i] += 1;

					RootWriter->AddHit(0,
									   1,
									   i,
									   offsets[i].first + hit1_rw.first,
									   offsets[i].second + hit1_rw.second,
									   hit1_rw.first,
									   hit1_rw.second,
									   smearedERing,
									   smearedEWedge,
									   localCoords.GetX(),
									   localCoords.GetY()
									   );

				}
			}


			/*//////////////////////////////////////////
			//  check if particle 2 (bu1) is detected //
			//////////////////////////////////////////*/
			smearedERing = 0.;
			smearedEWedge = 0.;

			double dtheta2 = straggler_par2_->Sample();//RIGHT NOW THIS CORRESPONDS TO 6LI GS IN LIF WORST CASE
			double dphi2 = straggler_par2_->SamplePhi();//RIGHT NOW THIS CORRESPONDS TO 6LI GS IN LIF WORST CASE

			//we now define the original kinematical trajectory as:
			Vec3 originalTrajectory2;
			originalTrajectory2.SetVectorSpherical(1, theta2*DEGRAD, phi2*DEGRAD);

			//define new basis vectors
			Vec3 etheta2, ephi2;
			etheta2.SetVectorCartesian(std::cos(theta2*DEGRAD)*std::cos(phi2*DEGRAD), std::cos(theta2*DEGRAD)*std::sin(phi2*DEGRAD), -std::sin(theta2*DEGRAD));
			ephi2.SetVectorCartesian(-std::sin(phi2*DEGRAD), std::cos(phi2*DEGRAD), 0.);

			//adjusted trajectory:
			Vec3 adjustedTrajectory2;
			adjustedTrajectory2 = std::cos(dtheta2*DEGRAD)*originalTrajectory2 + std::sin(dtheta2*DEGRAD)*(std::cos(dphi2*DEGRAD)*etheta2 + std::sin(dphi2*DEGRAD)*ephi2);
			adjustedTrajectory2 = adjustedTrajectory2.Unit();

			//now convert back
			double theta2_prime, phi2_prime;
			if(targetStraggle2){

				theta2_prime = adjustedTrajectory2.GetTheta()*RADDEG;//std::acos(adjustedTrajectory2.GetZ())*RADDEG;
				phi2_prime = adjustedTrajectory2.GetPhi()*RADDEG;//std::atan2(adjustedTrajectory2.GetY(), adjustedTrajectory2.GetX())*RADDEG;
				if(phi2_prime < 0) phi2_prime += 360.;
				RootPlotter->FillStraggleHistos(theta2, phi2, theta2_prime, phi2_prime, dtheta2, dphi2);

			} else {

				theta2_prime = theta2;
				phi2_prime = phi2;
				RootPlotter->FillStraggleHistos(theta2, phi2, theta2_prime, phi2_prime, 0, 0);

			}

			std::pair<int,int> hit2_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta2_prime*DEGRAD,phi2_prime*DEGRAD,reactionOrigin);

			if(hit2_rw.first != -1 && hit2_rw.second != -1 && !detected1){
				//std::cout << "here2" << std::endl;
				//apply energy loss to e1:
				double e2_aftertarget = targetLoss_par2_->ApplyEnergyLoss(e2, theta2_prime);

				//apply dead layer energy loss to e2_aftertarget
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta2_prime*DEGRAD,phi2_prime*DEGRAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e2_afterDeadLayer = deadLayerLoss_par2_->ApplyEnergyLoss(e2_aftertarget, trajectory, normal);

				//std::cout << "here2 again" << std::endl;

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit2_rw.first,e2_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit2_rw.second,e2_afterDeadLayer,smearedEWedge)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit2_rw.first,hit2_rw.second);
					ss << 200+i << "\t" << hit2_rw.first << "\t" << hit2_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected1 = true;
					hit1_ += 1;
					detectorHits_[i] += 1;

					RootWriter->AddHit(1,
									   2,
									   i,
									   offsets[i].first + hit2_rw.first,
									   offsets[i].second + hit2_rw.second,
									   hit2_rw.first,
									   hit2_rw.second,
									   smearedERing,
									   smearedEWedge,
									   localCoords.GetX(),
									   localCoords.GetY()
									   );
				}

				//std::cout << "here2 again again" << std::endl;
			}

			/*//////////////////////////////////////////
			//  check if particle 3 (bu2) is detected //
			//////////////////////////////////////////*/
			smearedERing = 0.;
			smearedEWedge = 0.;

			double dtheta3 = straggler_par3_->Sample();
			double dphi3 = straggler_par3_->SamplePhi();

			//define original kinematical trajectory
			Vec3 originalTrajectory3;
			originalTrajectory3.SetVectorSpherical(1, theta3*DEGRAD, phi3*DEGRAD);

			//define new basis vectors
			Vec3 etheta3, ephi3;
			etheta3.SetVectorCartesian(std::cos(theta3*DEGRAD)*std::cos(phi3*DEGRAD), std::cos(theta3*DEGRAD)*std::sin(phi3*DEGRAD), -std::sin(theta3*DEGRAD));
			ephi3.SetVectorCartesian(-std::sin(phi3*DEGRAD), std::cos(phi3*DEGRAD), 0.);

			//adjusted trajectory:
			Vec3 adjustedTrajectory3;
			adjustedTrajectory3 = std::cos(dtheta3*DEGRAD)*originalTrajectory3 + std::sin(dtheta3*DEGRAD)*(std::cos(dphi3*DEGRAD)*etheta3 + std::sin(dphi3*DEGRAD)*ephi3);
			adjustedTrajectory3 = adjustedTrajectory3.Unit();

			//now convert back
			double theta3_prime, phi3_prime;
			if(targetStraggle3){

				theta3_prime = adjustedTrajectory3.GetTheta()*RADDEG;
				phi3_prime = adjustedTrajectory3.GetPhi()*RADDEG;
				if(phi3_prime < 0) phi3_prime += 360.;
				RootPlotter->FillStraggleHistos(theta3, phi3, theta3_prime, phi3_prime, dtheta3, dphi3);

			} else {

				theta3_prime = theta3;
				phi3_prime = phi3;
				RootPlotter->FillStraggleHistos(theta3, phi3, theta3_prime, phi3_prime, dtheta3, dphi3);

			}

			std::pair<int,int> hit3_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta3_prime*DEGRAD,phi3_prime*DEGRAD,reactionOrigin);

			if(hit3_rw.first != -1 && hit3_rw.second != -1 && !detected2){
				//std::cout << "here3" << std::endl;
				//apply energy loss to e3:
				double e3_aftertarget = targetLoss_par3_->ApplyEnergyLoss(e3, theta3_prime);

				//apply dead layer energy loss to e3_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta3_prime*DEGRAD,phi3_prime*DEGRAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal * (1/normal.Mag());
				double e3_afterDeadLayer = deadLayerLoss_par3_->ApplyEnergyLoss(e3_aftertarget, trajectory, normal);

				//std::cout << "here3 again" << std::endl;

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit3_rw.first,e3_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit3_rw.second,e3_afterDeadLayer,smearedEWedge)){
					//std::cout << "here3 again again 1" << std::endl;
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit3_rw.first, hit3_rw.second);
					//std::cout << "here3 again again 2" << std::endl;
					ss << 300+i << "\t" << hit3_rw.first << "\t" << hit3_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					//std::cout << "here3 again again 3" << std::endl;
					detected2 = true;
					//std::cout << "here3 again again 4" << std::endl;
					hit2_ += 1;
					//std::cout << "here3 again again 5" << std::endl;
					detectorHits_[i] += 1;
					//std::cout << "here3 again again 6" << std::endl;

					RootWriter->AddHit(2,
									   3,
									   i,
									   offsets[i].first + hit3_rw.first,
									   offsets[i].second + hit3_rw.second,
									   hit3_rw.first,
									   hit3_rw.second,
									   smearedERing,
									   smearedEWedge,
									   localCoords.GetX(),
									   localCoords.GetY()
									   );
				}

				//std::cout << "here3 again again again" << std::endl;

			}

			/*//////////////////////////////////////////
			//  check if particle 4 (bu2) is detected //
			//////////////////////////////////////////*/
			smearedERing = 0.;
			smearedEWedge = 0.;

			double dtheta4 = straggler_par4_->Sample();
			double dphi4 = straggler_par4_->SamplePhi();

			//define original kinematical trajectory
			Vec3 originalTrajectory4;
			originalTrajectory4.SetVectorSpherical(1, theta4*DEGRAD, phi4*DEGRAD);

			//define new basis vectors
			Vec3 etheta4, ephi4;
			etheta4.SetVectorCartesian(std::cos(theta4*DEGRAD)*std::cos(phi4*DEGRAD), std::cos(theta4*DEGRAD)*std::sin(phi4*DEGRAD), -std::sin(theta4*DEGRAD));
			ephi4.SetVectorCartesian(-std::sin(phi4*DEGRAD), std::cos(phi4*DEGRAD), 0.);

			//adjusted trajectory:
			Vec3 adjustedTrajectory4;
			adjustedTrajectory4 = std::cos(dtheta4*DEGRAD)*originalTrajectory4 + std::sin(dtheta4*DEGRAD)*(std::cos(dphi4*DEGRAD)*etheta4 + std::sin(dphi4*DEGRAD)*ephi4);
			adjustedTrajectory4 = adjustedTrajectory4.Unit();

			//now convert back
			double theta4_prime, phi4_prime;
			if(targetStraggle4){

				theta4_prime = adjustedTrajectory4.GetTheta()*RADDEG;
				phi4_prime = adjustedTrajectory4.GetPhi()*RADDEG;
				if(phi4_prime < 0) phi4_prime += 360.;
				RootPlotter->FillStraggleHistos(theta4, phi4, theta4_prime, phi4_prime, dtheta4, dphi4);

			} else {

				theta4_prime = theta4;
				phi4_prime = phi4;
				RootPlotter->FillStraggleHistos(theta4, phi4, theta4_prime, phi4_prime, dtheta4, dphi4);

			}

			std::pair<int,int> hit4_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta4_prime*DEGRAD,phi4_prime*DEGRAD,reactionOrigin);

			if(hit4_rw.first != -1 && hit4_rw.second != -1 && !detected3){
				//std::cout << "here4" << std::endl;
				//apply energy loss to e4:
				double e4_aftertarget = targetLoss_par4_->ApplyEnergyLoss(e4, theta4_prime);

				//apply dead layer energy loss to e3_aftertarget
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta4_prime*DEGRAD,phi4_prime*DEGRAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal * (1/normal.Mag());
				double e4_afterDeadLayer = deadLayerLoss_par4_->ApplyEnergyLoss(e4_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit4_rw.first,e4_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit4_rw.second,e4_afterDeadLayer,smearedEWedge)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit4_rw.first,hit4_rw.second);
					ss << 400+i << "\t" << hit4_rw.first << "\t" << hit4_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected3 = true;
					hit3_ += 1;
					detectorHits_[i] += 1;

					RootWriter->AddHit(3,
									   4,
									   i,
									   offsets[i].first + hit4_rw.first,
									   offsets[i].second + hit4_rw.second,
									   hit4_rw.first,
									   hit4_rw.second,
									   smearedERing,
									   smearedEWedge,
									   localCoords.GetX(),
									   localCoords.GetY()
									   );
				}
			}

		//std::cout << "here end 1" << std::endl;

		}

		//std::cout << "here end 2" << std::endl;

		ss << eoev;// << std::endl;

		outfile << ss.str() << "\n";

		RootPlotter->ProcessTXTOutput(ss.str());

		RootWriter->FillEvent();

		// if(detectedEj) hitej_ +=  1;
		// if(detected1) hit1_ += 1;
		// if(detected2) hit2_ += 1;
		// if(detected3) hit3_ += 1;
		
		if(detectedEj && !detected1 && !detected2 && !detected3) hitOnlyEj_ += 1;
		if(!detectedEj && detected1 && !detected2 && !detected3) hitOnly1_ += 1;
		if(!detectedEj && !detected1 && detected2 && !detected3) hitOnly2_ += 1;
		if(!detectedEj && !detected1 && !detected2 && detected3) hitOnly3_ += 1;

		if(!detectedEj && detected1 && detected2 && !detected3) hitOnly12_ += 1;
		if(!detectedEj && detected1 && !detected2 && detected3) hitOnly13_ += 1;
		if(!detectedEj && !detected1 && detected2 && detected3) hitOnly23_ += 1;

		if(detected2&&detected3) hitBoth23_ += 1;


		if(!detectedEj && detected1 && detected2 && detected3) hitOnly123_ += 1;

		if((detectedEj && !detected1 && !detected2 && !detected3) || (!detectedEj && detected1 && !detected2 && !detected3) || (!detectedEj && !detected1 && detected2 && !detected3) || (!detectedEj && !detected1 && !detected2 && detected3)) onePartHits_ += 1;
		if((detectedEj && detected1 && !detected2 && !detected3) || (!detectedEj && detected1 && detected2 && !detected3) || (!detectedEj && !detected1 && detected2 && detected3) || (detectedEj && !detected1 && !detected2 && detected3)) twoPartHits_ += 1;
		if((detectedEj && detected1 && detected2 && !detected3) || (!detectedEj && detected1 && detected2 && detected3) || (detectedEj && !detected1 && detected2 && detected3) || (detectedEj && detected1 && !detected2 && detected3)) threePartHits_ += 1;
		if(detectedEj && detected1 && detected2 && detected3) fourPartHits_ += 1;

	}

	TFile *tempfile = new TFile("BeamSpotHisto_det3mc.root","RECREATE");
	//hBeamSpot->Write();
	tempfile->Close();
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
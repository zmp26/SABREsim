#include "det3mc.h"
#include <iostream>
#include "ConsoleColorizer.h"
#include "TH2.h"
#include "TFile.h"
#include "IMMMA_Tool_3.h"

static const int eoev = -1;//end of event value (eoev), printed between events in .det file (-1 will separate entries in mass text file)

const std::pair<int,int> det3mc::offsets[] = {
		{112,40},	//detector0 {ringOffset,wedgeOffset}
		{96,32},	//detector1 {ringOffset,wedgeOffset}
		{80,16},	//detector2 {ringOffset,wedgeOffset}
		{64,24},	//detector3 {ringOffset,wedgeOffset}
		{48,0}		//detector4 {ringOffset,wedgeOffset}
	};

det3mc::det3mc(std::vector<SABRE_Detector*>& SABRE_Array,
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
	: SABRE_Array_(SABRE_Array),
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
	  nevents_(0),
	  hit1_(0), hit3_(0), hit4_(0),
	  hitBoth34_(0), hitOnly3_(0), hitOnly4_(0),
	  onePartHits_(0), twoPartHits_(0), threePartHits_(0),
	  detectorHits_(SABRE_Array.size()),
	  beamspot_(beamspot)
	{

	}

// void det3mc::Run(std::ifstream& infile, std::ofstream& outfile, RootWriter* RootWriter, plot3mc* RootPlotter, bool targetStraggle1=true, bool targetStraggle2=true, bool targetStraggle3=true, bool targetStraggle4=true){
	
// 	double e1, theta1, phi1, e2, theta2, phi2, e3, theta3, phi3, e4, theta4, phi4;
	
// 	//TH2D *hBeamSpot = new TH2D("hBeamSpot","BeamSpot",200, -0.05, 0.05, 200, -0.05, 0.05);

// 	while(infile >> e1 >> theta1 >> phi1 >> e2 >> theta2 >> phi2 >> e3 >> theta3 >> phi3 >> e4 >> theta4 >> phi4){
		
// 		nevents_ += 1;
// 		bool detected1 = false, detected3 = false, detected4 = false;

// 		//get reaction origin based on beamspot
// 		Vec3 reactionOrigin = beamspot_->GeneratePoint();//same for whole event!
// 		//hBeamSpot->Fill(reactionOrigin.GetX(),reactionOrigin.GetY());
// 		//Vec3 reactionOrigin = {0.,0.,0.};

// 		std::ostringstream ss;
// 		ss << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\t" << e3 << "\t" << theta3 << "\t" << phi3 << "\t" << e4 << "\t" << theta4 << "\t" << phi4 << std::endl;

// 		if(nevents_%50000==0) ConsoleColorizer::PrintBlue("Processed " + std::to_string(nevents_) + " events...\n"); //std::cout << "Processed " << nevents_ << " events..." << std::endl;
		

// 		RootWriter->SetKinematics(0,e1,theta1,phi1);
// 		RootWriter->SetKinematics(1,e2,theta2,phi2);
// 		RootWriter->SetKinematics(2,e3,theta3,phi3);
// 		RootWriter->SetKinematics(3,e4,theta4,phi4);


// 		RootWriter->SetReactionOrigin(reactionOrigin.GetX(), reactionOrigin.GetY(), reactionOrigin.GetZ());


// 		for(size_t i=0; i<SABRE_Array_.size(); i++){
// 			/*///////////////////////////////////////////////
// 			//  check if particle 1 (ejectile) is detected //
// 			///////////////////////////////////////////////*/
// 			double smearedERing = 0., smearedEWedge = 0.;

// 			double dtheta1 = straggler_par1_->Sample();
// 			double dphi1 = straggler_par1_->SamplePhi();

// 			//original kinematic trajectory
// 			Vec3 originalTrajectory1;
// 			originalTrajectory1.SetVectorSpherical(1, theta1*DEGRAD, phi1*DEGRAD);
			
// 			//define new basis vectors
// 			Vec3 etheta1, ephi1;
// 			etheta1.SetVectorCartesian(std::cos(theta1*DEGRAD)*std::cos(phi1*DEGRAD), std::cos(theta1*DEGRAD)*std::sin(phi1*DEGRAD), -std::sin(theta1*DEGRAD));
// 			ephi1.SetVectorCartesian(-std::sin(phi1*DEGRAD), std::cos(phi1*DEGRAD), 0.);

// 			//adjusted trajectory:
// 			Vec3 adjustedTrajectory1;
// 			adjustedTrajectory1 = std::cos(dtheta1*DEGRAD)*originalTrajectory1 + std::sin(dtheta1*DEGRAD)*(std::cos(dphi1*DEGRAD)*etheta1 + std::sin(dphi1*DEGRAD)*ephi1);
// 			adjustedTrajectory1 = adjustedTrajectory1.Unit();

// 			//now convert back
// 			double theta1_prime, phi1_prime;
// 			if(targetStraggle1){

// 				theta1_prime = adjustedTrajectory1.GetTheta()*RADDEG;
// 				phi1_prime = adjustedTrajectory1.GetPhi()*RADDEG;
// 				if(phi1_prime < 0) phi1_prime += 360.;
// 				RootPlotter->FillStraggleHistos(theta1, phi1, theta1_prime, phi1_prime, dtheta1, dphi1);

// 			} else {

// 				theta1_prime = theta1;
// 				phi1_prime = phi1;
// 				RootPlotter->FillStraggleHistos(theta1, phi1, theta1_prime, phi1_prime, 0, 0);

// 			}

// 			//std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta1*DEGRAD,phi1*DEGRAD);
// 			std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta1_prime*DEGRAD,phi1_prime*DEGRAD,reactionOrigin);

// 			if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){

// 				//apply target energy loss to e1:
// 				double e1_aftertarget = targetLoss_par1_->ApplyEnergyLoss(e1, theta1_prime);

// 				//apply dead layer energy loss to e1_aftertarget:
// 				Vec3 trajectory;
// 				trajectory.SetVectorSpherical(1,theta1_prime*DEGRAD,phi1_prime*DEGRAD);
// 				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
// 				normal = normal*(1/normal.Mag());
// 				double e1_afterDeadLayer = deadLayerLoss_par1_->ApplyEnergyLoss(e1_aftertarget, trajectory, normal);

// 				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit1_rw.first,e1_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit1_rw.second,e1_afterDeadLayer,smearedEWedge)){
// 					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
// 					//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",100+i,hit1_rw.first,hit1_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
// 					ss << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
// 					detected1=true;
// 					hit1_+=1;
// 					detectorHits_[i] += 1;

// 					RootWriter->AddHit(0,
// 									   1,
// 									   i,
// 									   offsets[i].first + hit1_rw.first,
// 									   offsets[i].second + hit1_rw.second,
// 									   hit1_rw.first,
// 									   hit1_rw.second,
// 									   smearedERing,
// 									   smearedEWedge,
// 									   localCoords.GetX(),
// 									   localCoords.GetY()
// 									   );

// 				}
// 			}
// 			/*//////////////////////////////////////////
// 			//  check if particle 3 (bu1) is detected //
// 			//////////////////////////////////////////*/
// 			smearedERing = 0.;
// 			smearedEWedge = 0.;

// 			double dtheta3 = straggler_par3_->Sample();
// 			double dphi3 = straggler_par3_->SamplePhi();

// 			//define original kinematical trajectory
// 			Vec3 originalTrajectory3;
// 			originalTrajectory3.SetVectorSpherical(1, theta3*DEGRAD, phi3*DEGRAD);

// 			//define new basis vectors
// 			Vec3 etheta3, ephi3;
// 			etheta3.SetVectorCartesian(std::cos(theta3*DEGRAD)*std::cos(phi3*DEGRAD), std::cos(theta3*DEGRAD)*std::sin(phi3*DEGRAD), -std::sin(theta3*DEGRAD));
// 			ephi3.SetVectorCartesian(-std::sin(phi3*DEGRAD), std::cos(phi3*DEGRAD), 0.);

// 			//adjusted trajectory:
// 			Vec3 adjustedTrajectory3;
// 			adjustedTrajectory3 = std::cos(dtheta3*DEGRAD)*originalTrajectory3 + std::sin(dtheta3*DEGRAD)*(std::cos(dphi3*DEGRAD)*etheta3 + std::sin(dphi3*DEGRAD)*ephi3);
// 			adjustedTrajectory3 = adjustedTrajectory3.Unit();

// 			//now convert back
// 			double theta3_prime, phi3_prime;
// 			if(targetStraggle3){

// 				theta3_prime = adjustedTrajectory3.GetTheta()*RADDEG;
// 				phi3_prime = adjustedTrajectory3.GetPhi()*RADDEG;
// 				if(phi3_prime < 0) phi3_prime += 360.;
// 				RootPlotter->FillStraggleHistos(theta3, phi3, theta3_prime, phi3_prime, dtheta3, dphi3);

// 			} else {

// 				theta3_prime = theta3;
// 				phi3_prime = phi3;
// 				RootPlotter->FillStraggleHistos(theta3, phi3, theta3_prime, phi3_prime, dtheta3, dphi3);

// 			}

// 			//std::pair<int,int> hit3_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta3*DEGRAD,phi3*DEGRAD);
// 			std::pair<int,int> hit3_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta3_prime*DEGRAD,phi3_prime*DEGRAD,reactionOrigin);

// 			if(hit3_rw.first != -1 && hit3_rw.second != -1 && !detected3){

// 				//apply target energy loss to e3:
// 				double e3_aftertarget = targetLoss_par3_->ApplyEnergyLoss(e3, theta3_prime);

// 				//apply dead layer energy loss to e3_aftertarget:
// 				Vec3 trajectory;
// 				trajectory.SetVectorSpherical(1,theta3_prime*DEGRAD,phi3_prime*DEGRAD);
// 				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
// 				normal = normal*(1/normal.Mag());
// 				double e3_afterDeadLayer = deadLayerLoss_par3_->ApplyEnergyLoss(e3_aftertarget, trajectory, normal);

// 				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit3_rw.first,e3_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit3_rw.second,e3_afterDeadLayer,smearedEWedge)){
// 					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit3_rw.first,hit3_rw.second);
// 					//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",300+i,hit3_rw.first,hit3_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
// 					ss << 300+i << "\t" << hit3_rw.first << "\t" << hit3_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
// 					detected3 = true;
// 					hit3_ += 1;
// 					detectorHits_[i] += 1;

// 					RootWriter->AddHit(2,
// 									   3,
// 									   i,
// 									   offsets[i].first + hit3_rw.first,
// 									   offsets[i].second + hit3_rw.second,
// 									   hit3_rw.first,
// 									   hit3_rw.second,
// 									   smearedERing,
// 									   smearedEWedge,
// 									   localCoords.GetX(),
// 									   localCoords.GetY()
// 									   );

// 				}
// 			}


// 			/*//////////////////////////////////////////
// 			//  check if particle 4 (bu2) is detected //
// 			//////////////////////////////////////////*/
// 			smearedERing = 0.;
// 			smearedEWedge = 0.;

// 			double dtheta4 = straggler_par4_->Sample();
// 			double dphi4 = straggler_par4_->SamplePhi();

// 			//define original kinematical trajectory
// 			Vec3 originalTrajectory4;
// 			originalTrajectory4.SetVectorSpherical(1, theta4*DEGRAD, phi4*DEGRAD);

// 			//define new basis vectors
// 			Vec3 etheta4, ephi4;
// 			etheta4.SetVectorCartesian(std::cos(theta4*DEGRAD)*std::cos(phi4*DEGRAD), std::cos(theta4*DEGRAD)*std::sin(phi4*DEGRAD), -std::sin(theta4*DEGRAD));
// 			ephi4.SetVectorCartesian(-std::sin(phi4*DEGRAD), std::cos(phi4*DEGRAD), 0.);

// 			//adjusted trajectory:
// 			Vec3 adjustedTrajectory4;
// 			adjustedTrajectory4 = std::cos(dtheta4*DEGRAD)*originalTrajectory4 + std::sin(dtheta4*DEGRAD)*(std::cos(dphi4*DEGRAD)*etheta4 + std::sin(dphi4*DEGRAD)*ephi4);
// 			adjustedTrajectory4 = adjustedTrajectory4.Unit();

// 			//now convert back
// 			double theta4_prime, phi4_prime;
// 			if(targetStraggle4){

// 				theta4_prime = adjustedTrajectory4.GetTheta()*RADDEG;
// 				phi4_prime = adjustedTrajectory4.GetPhi()*RADDEG;
// 				if(phi4_prime < 0) phi4_prime += 360.;
// 				RootPlotter->FillStraggleHistos(theta4, phi4, theta4_prime, phi4_prime, dtheta4, dphi4);

// 			} else {

// 				theta4_prime = theta4;
// 				phi4_prime = phi4;
// 				RootPlotter->FillStraggleHistos(theta4, phi4, theta4_prime, phi4_prime, dtheta4, dphi4);

// 			}

// 			//std::pair<int,int> hit4_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta4*DEGRAD,phi4*DEGRAD);
// 			std::pair<int,int> hit4_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta4_prime*DEGRAD,phi4_prime*DEGRAD,reactionOrigin);

// 			if(hit4_rw.first != -1 && hit4_rw.second != -1 && !detected4){

// 				//apply target energy loss to e4:
// 				double e4_aftertarget = targetLoss_par4_->ApplyEnergyLoss(e4, theta4_prime);

// 				//apply dead layer energy loss to e4_aftertarget:
// 				Vec3 trajectory;
// 				trajectory.SetVectorSpherical(1,theta4_prime*DEGRAD,phi4_prime*DEGRAD);
// 				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
// 				normal = normal*(1/normal.Mag());
// 				double e4_afterDeadLayer = deadLayerLoss_par4_->ApplyEnergyLoss(e4_aftertarget, trajectory, normal);

// 				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit4_rw.first,e4_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit4_rw.second,e4_afterDeadLayer,smearedEWedge)){
// 					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit4_rw.first,hit4_rw.second);
// 					//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",400+i,hit4_rw.first,hit4_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
// 					ss << 400+i << "\t" << hit4_rw.first << "\t" << hit4_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
// 					detected4 = true;
// 					hit4_ += 1;
// 					detectorHits_[i] += 1;

// 					RootWriter->AddHit(3,
// 									   4,
// 									   i,
// 									   offsets[i].first + hit4_rw.first,
// 									   offsets[i].second + hit4_rw.second,
// 									   hit4_rw.first,
// 									   hit4_rw.second,
// 									   smearedERing,
// 									   smearedEWedge,
// 									   localCoords.GetX(),
// 									   localCoords.GetY()
// 									   );
// 				}
// 			}

// 		}

// 		if(detected3 && detected4){
// 			hitBoth34_ += 1;
// 		}
// 		if(detected3&&!detected4){
// 			hitOnly3_ += 1;
// 		}
// 		if(detected4&&!detected3){
// 			hitOnly4_ += 1;
// 		}

// 		if((detected1&&!detected3&&!detected4) || (!detected1&&detected3&&!detected4) || (!detected1&&!detected3&&detected4)){
// 			onePartHits_ += 1;
// 		} else if((detected1&&detected3&&!detected4) || (!detected1&&detected3&&detected4) || (detected1&&!detected3&&detected4)){
// 			twoPartHits_ += 1;
// 		} else if((detected1&&detected3&&detected4)){
// 			threePartHits_ += 1;
// 		}

// 		ss << eoev;

// 		outfile << ss.str() << "\n";

// 		RootPlotter->ProcessTXTOutput(ss.str());

// 		RootWriter->FillEvent();
// 	}

// 	TFile *tempfile = new TFile("BeamSpotHisto_det3mc.root","RECREATE");
// 	//hBeamSpot->Write();
// 	tempfile->Close();
// }

void det3mc::Run(std::ifstream& infile, std::ofstream& outfile, RootWriter* RootWriter, plot3mc* RootPlotter, bool targetStraggle1=true, bool targetStraggle2=true, bool targetStraggle3=true, bool targetStraggle4=true){
	
	double e1, theta1, phi1, e2, theta2, phi2, e3, theta3, phi3, e4, theta4, phi4;
	
	//TH2D *hBeamSpot = new TH2D("hBeamSpot","BeamSpot",200, -0.05, 0.05, 200, -0.05, 0.05);


	// ----- Initialize IMMMA_Tool_3 Below -----
	// -----	This should be updated     -----
	// -----	before running, and will   -----
	// -----	eventually come from       -----
	// -----	SimConfig, but that will   -----
	// -----	take some more time        -----
	IMMMA_Tool_3 immmaTool;
	immmaTool.SetBeamNucleus(0,"",0);
	immmaTool.SetTargetNucleus(0,"",0);
	immmaTool.SetEjectileNucleus(0,"",0);
	immmaTool.SetRecoilNucleus(0,"",0);
	immmaTool.SetBreakup1Nucleus(0,"",0);
	immmaTool.SetBreakup2Nucleus(0,"",0);
	immmaTool.SetBeamEnergyMeV(0);
	immmaTool.SetRecoilExEMeV(0);


	while(infile >> e1 >> theta1 >> phi1 >> e2 >> theta2 >> phi2 >> e3 >> theta3 >> phi3 >> e4 >> theta4 >> phi4){
		
		nevents_ += 1;

		//get reaction origin based on beamspot
		Vec3 reactionOrigin = beamspot_->GeneratePoint();//same for whole event!
		//hBeamSpot->Fill(reactionOrigin.GetX(),reactionOrigin.GetY());
		//Vec3 reactionOrigin = {0.,0.,0.};

		std::ostringstream ss;
		ss << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\t" << e3 << "\t" << theta3 << "\t" << phi3 << "\t" << e4 << "\t" << theta4 << "\t" << phi4 << std::endl;

		if(nevents_%50000==0) ConsoleColorizer::PrintBlue("Processed " + std::to_string(nevents_) + " events...\n"); //std::cout << "Processed " << nevents_ << " events..." << std::endl;
		

		RootWriter->SetKinematics(0,e1,theta1,phi1);
		RootWriter->SetKinematics(1,e2,theta2,phi2);
		RootWriter->SetKinematics(2,e3,theta3,phi3);
		RootWriter->SetKinematics(3,e4,theta4,phi4);


		RootWriter->SetReactionOrigin(reactionOrigin.GetX(), reactionOrigin.GetY(), reactionOrigin.GetZ());


		double smearedERing1 = 0., smearedERing3 = 0., smearedERing4 = 0.;
		double smearedEWedge1 = 0., smearedEWedge3 = 0., smearedEWedge4 = 0.;
		double theta1_meas = 0., theta3_meas = 0., theta4_meas = 0.;
		double phi1_meas = 0., phi3_meas = 0., phi4_meas = 0.;
		bool detected1 = false, detected3 = false, detected4 = false;

		for(size_t i=0; i<SABRE_Array_.size(); i++){
			/*///////////////////////////////////////////////
			//  check if particle 1 (ejectile) is detected //
			///////////////////////////////////////////////*/

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

			//std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta1*DEGRAD,phi1*DEGRAD);
			std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta1_prime*DEGRAD,phi1_prime*DEGRAD,reactionOrigin);

			if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){

				//apply target energy loss to e1:
				double e1_aftertarget = targetLoss_par1_->ApplyEnergyLoss(e1, theta1_prime);

				//apply dead layer energy loss to e1_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta1_prime*DEGRAD,phi1_prime*DEGRAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e1_afterDeadLayer = deadLayerLoss_par1_->ApplyEnergyLoss(e1_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit1_rw.first,e1_afterDeadLayer,smearedERing1) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit1_rw.second,e1_afterDeadLayer,smearedEWedge1)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
					//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",100+i,hit1_rw.first,hit1_rw.second,smearedERing1,smearedEWedge1,localCoords.GetX(),localCoords.GetY()) << endl;
					ss << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing1 << "\t" << smearedEWedge1 << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected1=true;
					hit1_+=1;
					detectorHits_[i] += 1;

					RootWriter->AddHit(0,
									   1,
									   i,
									   offsets[i].first + hit1_rw.first,
									   offsets[i].second + hit1_rw.second,
									   hit1_rw.first,
									   hit1_rw.second,
									   smearedERing1,
									   smearedEWedge1,
									   localCoords.GetX(),
									   localCoords.GetY()
									   );

					//for now, update theta1_meas to be exactly theta1_prime -> eventually, add SABREAngleMap class and add functionality here and in plot2/3/4mc which currently have their own readAngleMaps
					//likewise for phi1_meas
					theta1_meas = theta1_prime;
					phi1_meas = phi1_prime;

				}
			}
			/*//////////////////////////////////////////
			//  check if particle 3 (bu1) is detected //
			//////////////////////////////////////////*/

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

			//std::pair<int,int> hit3_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta3*DEGRAD,phi3*DEGRAD);
			std::pair<int,int> hit3_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta3_prime*DEGRAD,phi3_prime*DEGRAD,reactionOrigin);

			if(hit3_rw.first != -1 && hit3_rw.second != -1 && !detected3){

				//apply target energy loss to e3:
				double e3_aftertarget = targetLoss_par3_->ApplyEnergyLoss(e3, theta3_prime);

				//apply dead layer energy loss to e3_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta3_prime*DEGRAD,phi3_prime*DEGRAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e3_afterDeadLayer = deadLayerLoss_par3_->ApplyEnergyLoss(e3_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit3_rw.first,e3_afterDeadLayer,smearedERing3) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit3_rw.second,e3_afterDeadLayer,smearedEWedge3)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit3_rw.first,hit3_rw.second);
					//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",300+i,hit3_rw.first,hit3_rw.second,smearedERing3,smearedEWedge3,localCoords.GetX(),localCoords.GetY()) << endl;
					ss << 300+i << "\t" << hit3_rw.first << "\t" << hit3_rw.second << "\t" << smearedERing3 << "\t" << smearedEWedge3 << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected3 = true;
					hit3_ += 1;
					detectorHits_[i] += 1;

					RootWriter->AddHit(2,
									   3,
									   i,
									   offsets[i].first + hit3_rw.first,
									   offsets[i].second + hit3_rw.second,
									   hit3_rw.first,
									   hit3_rw.second,
									   smearedERing3,
									   smearedEWedge3,
									   localCoords.GetX(),
									   localCoords.GetY()
									   );

					//for now, update theta3_meas to be exactly theta3_prime -> eventually, add SABREAngleMap class and add functionality here and in plot2/3/4mc which currently have their own readAngleMaps
					//likewise for phi3_meas
					theta3_meas = theta3_prime;
					phi3_meas = phi3_prime;

				}
			}


			/*//////////////////////////////////////////
			//  check if particle 4 (bu2) is detected //
			//////////////////////////////////////////*/

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

			//std::pair<int,int> hit4_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta4*DEGRAD,phi4*DEGRAD);
			std::pair<int,int> hit4_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta4_prime*DEGRAD,phi4_prime*DEGRAD,reactionOrigin);

			if(hit4_rw.first != -1 && hit4_rw.second != -1 && !detected4){

				//apply target energy loss to e4:
				double e4_aftertarget = targetLoss_par4_->ApplyEnergyLoss(e4, theta4_prime);

				//apply dead layer energy loss to e4_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta4_prime*DEGRAD,phi4_prime*DEGRAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e4_afterDeadLayer = deadLayerLoss_par4_->ApplyEnergyLoss(e4_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit4_rw.first,e4_afterDeadLayer,smearedERing4) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit4_rw.second,e4_afterDeadLayer,smearedEWedge4)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit4_rw.first,hit4_rw.second);
					//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",400+i,hit4_rw.first,hit4_rw.second,smearedERing4,smearedEWedge4,localCoords.GetX(),localCoords.GetY()) << endl;
					ss << 400+i << "\t" << hit4_rw.first << "\t" << hit4_rw.second << "\t" << smearedERing4 << "\t" << smearedEWedge4 << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected4 = true;
					hit4_ += 1;
					detectorHits_[i] += 1;

					RootWriter->AddHit(3,//hit index
									   4,//particle ID
									   i,//detector id
									   offsets[i].first + hit4_rw.first,//ring channel
									   offsets[i].second + hit4_rw.second,//wedge channel
									   hit4_rw.first,//local ring
									   hit4_rw.second,//local wedge
									   smearedERing4,//ring energy
									   smearedEWedge4,//wedge energy
									   localCoords.GetX(),//local x
									   localCoords.GetY()//local y
									   );

					//for now, update theta2_meas to be exactly theta2_prime -> eventually, add SABREAngleMap class and add functionality here and in plot2/3/4mc which currently have their own readAngleMaps
					//likewise for phi2_meas
					theta4_meas = theta4_prime;
					phi4_meas = phi4_prime;
				}
			}

		}

		if(detected3 && detected4){
			//both break up fragments detected
			hitBoth34_ += 1;

			//IMM Analysis:
			auto immmaResults = immmaTool.AnalyzeEventIMM(smearedERing1, theta1_meas, phi1_meas,
														  smearedERing3, theta3_meas, phi3_meas,
														  smearedERing4, theta4_meas, phi4_meas);

			const CaseResult& resultA = immmaResults.first;
			const CaseResult& resultB = immmaResults.second;

			//pass resultA, resultB into RootWriter and RootPlotter here eventually, but must add that functionality first!

		} else if(detected3 || detected4){
			//only one breakup fragment detected
			double detectedE_meas = detected3 ? smearedERing3 : smearedERing4;
			double detectedTheta_meas = detected3 ? theta3_meas : theta4_meas;
			double detectedPhi_meas = detected3 ? phi3_meas : phi4_meas;

			auto immmaResults = immmaTool.AnalyzeEventMMM(smearedERing1, theta1_meas, phi1_meas,
														  detectedE_meas, detectedTheta_meas, detectedPhi_meas);
			const CaseResult& resultA = immmaResults.first;
			const CaseResult& resultB = immmaResults.second;

			//pass resultA, resultB into RootWriter and RootPlotter here eventually, but must add that functionality first!

		}

		if((detected1&&!detected3&&!detected4) || (!detected1&&detected3&&!detected4) || (!detected1&&!detected3&&detected4)){
			onePartHits_ += 1;
		} else if((detected1&&detected3&&!detected4) || (!detected1&&detected3&&detected4) || (detected1&&!detected3&&detected4)){
			twoPartHits_ += 1;
		} else if((detected1&&detected3&&detected4)){
			threePartHits_ += 1;
		}

		ss << eoev;

		outfile << ss.str() << "\n";

		RootPlotter->ProcessTXTOutput(ss.str());

		RootWriter->FillEvent();
	}

	TFile *tempfile = new TFile("BeamSpotHisto_det3mc.root","RECREATE");
	//hBeamSpot->Write();
	tempfile->Close();
}

long det3mc::GetNumEvents() const {
	return nevents_;
}

long det3mc::GetHit1() const {
	return hit1_;
}

long det3mc::GetHit3() const {
	return hit3_;
}

long det3mc::GetHit4() const {
	return hit4_;
}

long det3mc::GetHitBoth34() const {
	return hitBoth34_;
}

long det3mc::GetHitOnly3() const {
	return hitOnly3_;
}

long det3mc::GetHitOnly4() const {
	return hitOnly4_;
}

long det3mc::GetOnePartHits() const {
	return onePartHits_;
}

long det3mc::GetTwoPartHits() const {
	return twoPartHits_;
}

long det3mc::GetThreePartHits() const {
	return threePartHits_;
}

const std::vector<long>& det3mc::GetDetectorHits() const {
	return detectorHits_;
}
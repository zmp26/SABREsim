#include "det2mc.h"
#include <iostream>
#include "ConsoleColorizer.h"
#include "TH2.h"
#include "TFile.h"

static const int eoev = -1;//end of event value (eoev), printed between events in .det file (-1 will separate entries in mass text file)

const std::pair<int,int> det2mc::offsets[] = {
		{112,40},	//detector0 {ringOffset,wedgeOffset}
		{96,32},	//detector1 {ringOffset,wedgeOffset}
		{80,16},	//detector2 {ringOffset,wedgeOffset}
		{64,24},	//detector3 {ringOffset,wedgeOffset}
		{48,0}		//detector4 {ringOffset,wedgeOffset}
	};

det2mc::det2mc(std::vector<SABRE_Detector*>& SABRE_Array,
			   std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels,
			   TargetEnergyLoss* targetLoss_par1,
			   TargetEnergyLoss* targetLoss_par2,
			   SABRE_DeadLayerModel* deadLayerLoss_par1,
			   SABRE_DeadLayerModel* deadLayerLoss_par2,
			   Beamspot* beamspot)
	: SABRE_Array_(SABRE_Array),
	  SABREARRAY_EnergyResolutionModels_(SABREARRAY_EnergyResolutionModels),
	  targetLoss_par1_(targetLoss_par1),
	  targetLoss_par2_(targetLoss_par2),
	  deadLayerLoss_par1_(deadLayerLoss_par1),
	  deadLayerLoss_par2_(deadLayerLoss_par2),
	  nevents_(0),
	  hit1_(0), hit2_(0), hitBoth_(0),
	  hit1Only_(0), hit2Only_(0),
	  detectorHits_(SABRE_Array.size(),0),
	  beamspot_(beamspot)
	  {

	  }


void det2mc::Run(std::ifstream& infile, std::ofstream& outfile, RootWriter* RootWriter, plot2mc* RootPlotter){
	double e1, theta1, phi1, thetacm, e2, theta2, phi2;

	//TH1D *hDeadLayerELoss = new TH1D("hDeadLayerELoss","hDeadLayerELoss;Energy (keV)", 30, 25, 28);
	TH2D *hBeamSpot = new TH2D("hBeamSpot","BeamSpot",200, -0.05, 0.05, 200, -0.05, 0.05);

	while(infile >> e1 >> theta1 >> phi1 >> thetacm >> e2 >> theta2 >> phi2){
		nevents_ += 1;
		bool detected1 = false;
		bool detected2 = false;

		//get reaction origin based on beamspot
		Vec3 reactionOrigin = beamspot_->GeneratePoint();//same for whole event!
		hBeamSpot->Fill(reactionOrigin.GetX(),reactionOrigin.GetY());
		//Vec3 reactionOrigin = {0.,0.,0.};

		std::ostringstream ss;
		//outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << thetacm << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\n";
		ss << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << thetacm << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\n";
		if(nevents_ % 50000 == 0) ConsoleColorizer::PrintBlue("Processed " + std::to_string(nevents_) + " events...\n");
		//outfile << ss.str();

		RootWriter->SetKinematics(0, e1, theta1, phi1);
		RootWriter->SetKinematics(1, e2, theta2, phi2);
		RootWriter->SetKinematics(2, -666., -666., -666.);
		RootWriter->SetKinematics(3, -666., thetacm, -666.);

		RootWriter->SetReactionOrigin(reactionOrigin.GetX(), reactionOrigin.GetY(), reactionOrigin.GetZ());

		for(size_t i=0; i<SABRE_Array_.size(); i++){
				/*///////////////////////////////////////////////
				//  check if particle 1 (ejectile) is detected //
				///////////////////////////////////////////////*/

				//prep variables to hold "smeared" ring and wedge energy after applying resolution
				double smearedERing = 0., smearedEWedge = 0.;

				//get <ring,wedge> pair based on theta,phi
				//std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD);
				std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta1*DEG2RAD, phi1*DEG2RAD, reactionOrigin);
				//pair<int,int> hit1_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD,{0.000001,0.000001,0.000001});
				if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){

					//apply target energy loss to e1:
					double e1_aftertarget = targetLoss_par1_->ApplyEnergyLoss(e1, theta1);

					//apply dead layer energy loss to e1_aftertarget:
					Vec3 trajectory;
					trajectory.SetVectorSpherical(1,theta1*DEG2RAD,phi1*DEG2RAD);
					Vec3 normal = SABRE_Array_[i]->GetNormTilted();
					normal = normal*(1/normal.Mag());
					double e1_afterDeadLayer = deadLayerLoss_par1_->ApplyEnergyLoss(e1_aftertarget, trajectory, normal);
					//hDeadLayerELoss->Fill(abs(e1_afterDeadLayer - e1_aftertarget)*1000.);

					// if(nevents%10000 == 0){
					// 	std::cout << "hit1 target_energy_loss = " << abs(e1-e1_aftertarget) << " MeV" << std::endl;
					// }

					if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit1_rw.first,e1_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit1_rw.second,e1_afterDeadLayer,smearedEWedge)){
						Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
						//outfile << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << "\n";
						ss << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << "\n";
						detected1 = true;
						hit1_+=1;
						detectorHits_[i] += 1;

						RootWriter->AddHit(0,					//hit index, 0 for first particle
										   1,					//particleID
										   i,					//detector id
										   offsets[i].first + hit1_rw.first,		//ring channel
										   offsets[i].second + hit1_rw.second,		//wedge channel
										   hit1_rw.first,		//local ring
										   hit1_rw.second,		//local wedge
										   smearedERing,
										   smearedEWedge,
										   localCoords.GetX(),
										   localCoords.GetY()
										   );



					}
				}


				/*/////////////////////////////////////////////
				//  check if particle 2 (recoil) is detected //
				/////////////////////////////////////////////*/

				smearedERing = 0.;
				smearedEWedge = 0.;
				//std::pair<int,int> hit2_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta2*DEG2RAD,phi2*DEG2RAD);
				std::pair<int,int> hit2_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta2*DEG2RAD,phi2*DEG2RAD,reactionOrigin);
				//pair<int,int> hit2_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta2*DEG2RAD,phi2*DEG2RAD,{0.000001,0.000001,0.000001});
				if(hit2_rw.first != -1 && hit2_rw.second != -1 && !detected2){

					//apply target energy loss to e1:
					double e2_aftertarget = targetLoss_par2_->ApplyEnergyLoss(e2, theta2);

					//apply dead layer energy loss to e1_aftertarget:
					Vec3 trajectory;
					trajectory.SetVectorSpherical(1,theta2*DEG2RAD,phi2*DEG2RAD);
					Vec3 normal = SABRE_Array_[i]->GetNormTilted();
					normal = normal*(1/normal.Mag());
					double e2_afterDeadLayer = deadLayerLoss_par2_->ApplyEnergyLoss(e2_aftertarget, trajectory, normal);
					//hDeadLayerELoss->Fill(std::fabs(e2_afterDeadLayer - e2_aftertarget)*1000.);

					//if(nevents%10000 == 0) std::cout << "hit2 target_energy_loss = " << e2-e2_aftertarget << " MeV" << std::endl;
					//hEnergyLosses->Fill(abs(e2-e2_aftertarget)*1000.);//for target energy loss
					//hEnergyLosses->Fill((e2_afterDeadLayer-e2_aftertarget)*1000.);//for dead layer energy loss

					if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit2_rw.first,e2_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit2_rw.second,e2_afterDeadLayer,smearedEWedge)){
						Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit2_rw.first,hit2_rw.second);
						//outfile << 200+i << "\t" << hit2_rw.first << "\t" << hit2_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << "\n";
						ss << 200+i << "\t" << hit2_rw.first << "\t" << hit2_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << "\n";
						detected2 = true;
						hit2_+=1;
						detectorHits_[i] += 1;

						RootWriter->AddHit(1,					//hit index, 0 for first particle
										   2,					//particleID
										   i,					//detector id
										   offsets[i].first + hit2_rw.first,		//ring channel
										   offsets[i].second + hit2_rw.second,		//wedge channel
										   hit2_rw.first,		//local ring
										   hit2_rw.second,		//local wedge
										   smearedERing,
										   smearedEWedge,
										   localCoords.GetX(),
										   localCoords.GetY()
										   );
					}
				}
		}

		if(detected1&&detected2) hitBoth_ += 1;

		if(detected1&&!detected2) hit1Only_ += 1;

		if(!detected1&&detected2) hit2Only_ += 1;

		ss << eoev;

		outfile << ss.str() << "\n";

		//std::cout << "ss.str() = " << ss.str() << "\n";
		RootPlotter->ProcessTXTOutput(ss.str());

		RootWriter->FillEvent();
	}

	// TFile *tempfile = new TFile("EnergyLossHisto.root","RECREATE");
	// //tempfile->cd();
	// hDeadLayerELoss->Write();
	// tempfile->Close();
	//std::cout << "here" << std::endl;

	// RootWriter->Set_detmc(2);
	// RootWriter->SetReaction();

	TFile *tempfile = new TFile("BeamSpotHisto_det2mc.root","RECREATE");
	hBeamSpot->Write();
	tempfile->Close();

}

long det2mc::GetNumEvents() const {
	return nevents_;
}

long det2mc::GetHit1() const {
	return hit1_;
}

long det2mc::GetHit2() const {
	return hit2_;
}

long det2mc::GetHitBoth() const {
	return hitBoth_;
}

long det2mc::GetHit1Only() const {
	return hit1Only_;
}

long det2mc::GetHit2Only() const {
	return hit2Only_;
}

const std::vector<long>& det2mc::GetDetectorHits() const {
	return detectorHits_;
}
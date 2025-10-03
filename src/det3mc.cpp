#include "det3mc.h"
#include <iostream>
#include "ConsoleColorizer.h"

static const int eoev = -1;//end of event value (eoev), printed between events in .det file (-1 will separate entries in mass text file)

det3mc::det3mc(std::vector<SABRE_Detector*>& SABRE_Array,
			   std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels,
			   TargetEnergyLoss* targetLoss,
			   SABRE_DeadLayerModel* deadLayerLoss)
	: SABRE_Array_(SABRE_Array),
	  SABREARRAY_EnergyResolutionModels_(SABREARRAY_EnergyResolutionModels),
	  targetLoss_(targetLoss),
	  deadLayerLoss_(deadLayerLoss),
	  nevents_(0),
	  hit1_(0), hit3_(0), hit4_(0),
	  hitBoth34_(0), hitOnly3_(0), hitOnly4_(0),
	  onePartHits_(0), twoPartHits_(0), threePartHits_(0),
	  detectorHits_(SABRE_Array.size())
	{

	}

void det3mc::Run(std::ifstream& infile, std::ofstream& outfile){
	double e1, theta1, phi1, e2, theta2, phi2, e3, theta3, phi3, e4, theta4, phi4;
	while(infile >> e1 >> theta1 >> phi1 >> e2 >> theta2 >> phi2 >> e3 >> theta3 >> phi3 >> e4 >> theta4 >> phi4){
		
		nevents_ += 1;
		bool detected1 = false, detected3 = false, detected4 = false;

		outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\t" << e3 << "\t" << theta3 << "\t" << phi3 << "\t" << e4 << "\t" << theta4 << "\t" << phi4 << std::endl;

		if(nevents_%50000==0) std::cout << "Processed " << nevents_ << " events..." << std::endl;

		for(size_t i=0; i<SABRE_Array_.size(); i++){
			/*///////////////////////////////////////////////
			//  check if particle 1 (ejectile) is detected //
			///////////////////////////////////////////////*/
			double smearedERing = 0., smearedEWedge = 0.;

			std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD);
			//std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD,{0., 0., 0.});

			if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){

				//apply target energy loss to e1:
				double e1_aftertarget = targetLoss_->ApplyEnergyLoss(e1, theta1);

				//apply dead layer energy loss to e1_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta1*DEG2RAD,phi1*DEG2RAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e1_afterDeadLayer = deadLayerLoss_->ApplyEnergyLoss(e1_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit1_rw.first,e1_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit1_rw.second,e1_afterDeadLayer,smearedEWedge)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
					//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",100+i,hit1_rw.first,hit1_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
					outfile << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected1=true;
					hit1_+=1;
					detectorHits_[i] += 1;
				}
			}
			/*//////////////////////////////////////////
			//  check if particle 3 (bu1) is detected //
			//////////////////////////////////////////*/
			smearedERing = 0.;
			smearedEWedge = 0.;

			std::pair<int,int> hit3_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta3*DEG2RAD,phi3*DEG2RAD);
			//std::pair<int,int> hit3_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta3*DEG2RAD,phi3*DEG2RAD,{0., 0., 0.});

			if(hit3_rw.first != -1 && hit3_rw.second != -1 && !detected3){

				//apply target energy loss to e3:
				double e3_aftertarget = targetLoss_->ApplyEnergyLoss(e3, theta3);

				//apply dead layer energy loss to e3_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta3*DEG2RAD,phi3*DEG2RAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e3_afterDeadLayer = deadLayerLoss_->ApplyEnergyLoss(e3_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit3_rw.first,e3_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit3_rw.second,e3_afterDeadLayer,smearedEWedge)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit3_rw.first,hit3_rw.second);
					//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",300+i,hit3_rw.first,hit3_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
					outfile << 300+i << "\t" << hit3_rw.first << "\t" << hit3_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected3 = true;
					hit3_ += 1;
					detectorHits_[i] += 1;
				}
			}


			/*//////////////////////////////////////////
			//  check if particle 4 (bu2) is detected //
			//////////////////////////////////////////*/
			smearedERing = 0.;
			smearedEWedge = 0.;

			std::pair<int,int> hit4_rw = SABRE_Array_[i]->GetTrajectoryRingWedge(theta4*DEG2RAD,phi4*DEG2RAD);
			//std::pair<int,int> hit4_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta4*DEG2RAD,phi4*DEG2RAD,{0., 0., 0.});

			if(hit4_rw.first != -1 && hit4_rw.second != -1 && !detected4){

				//apply target energy loss to e4:
				double e4_aftertarget = targetLoss_->ApplyEnergyLoss(e4, theta4);

				//apply dead layer energy loss to e4_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta4*DEG2RAD,phi4*DEG2RAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e4_afterDeadLayer = deadLayerLoss_->ApplyEnergyLoss(e4_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit4_rw.first,e4_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit4_rw.second,e4_afterDeadLayer,smearedEWedge)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit4_rw.first,hit4_rw.second);
					//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",400+i,hit4_rw.first,hit4_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
					outfile << 400+i << "\t" << hit4_rw.first << "\t" << hit4_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected4 = true;
					hit4_ += 1;
					detectorHits_[i] += 1;
				}
			}

		}

		outfile << eoev << std::endl;

		if(detected3 && detected4){
			hitBoth34_ += 1;
		}
		if(detected3&&!detected4){
			hitOnly3_ += 1;
		}
		if(detected4&&!detected3){
			hitOnly4_ += 1;
		}

		if((detected1&&!detected3&&!detected4) || (!detected1&&detected3&&!detected4) || (!detected1&&!detected3&&detected4)){
			onePartHits_ += 1;
		} else if((detected1&&detected3&&!detected4) || (!detected1&&detected3&&detected4) || (detected1&&!detected3&&detected4)){
			twoPartHits_ += 1;
		} else if((detected1&&detected3&&detected4)){
			threePartHits_ += 1;
		}
	}
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
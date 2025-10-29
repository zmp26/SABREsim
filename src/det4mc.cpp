#include "det4mc.h"
#include <iostream>
#include "ConsoleColorizer.h"

static const int eoev = -1;//end of event value (eoev), printed between events in .det file (-1 will separate entries in mass text file)

det4mc::det4mc(std::vector<SABRE_Detector*>& SABRE_Array,
			   std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels,
			   TargetEnergyLoss* targetLoss_par1, TargetEnergyLoss* targetLoss_par2, TargetEnergyLoss* targetLoss_par3, TargetEnergyLoss* targetLoss_par4,
			   SABRE_DeadLayerModel* deadLayerLoss_par1, SABRE_DeadLayerModel* deadLayerLoss_par2, SABRE_DeadLayerModel* deadLayerLoss_par3, SABRE_DeadLayerModel* deadLayerLoss_par4,
			   Beamspot* beamspot)
	: SABRE_Array_(SABRE_Array),
	  SABREARRAY_EnergyResolutionModels_(SABREARRAY_EnergyResolutionModels),
	  targetLoss_par1_(targetLoss_par1), targetLoss_par2_(targetLoss_par2), targetLoss_par3_(targetLoss_par3), targetLoss_par4_(targetLoss_par4),
	  deadLayerLoss_par1_(deadLayerLoss_par1), deadLayerLoss_par2_(deadLayerLoss_par2), deadLayerLoss_par3_(deadLayerLoss_par3), deadLayerLoss_par4_(deadLayerLoss_par4),
	  nevents_(0),
	  hitej_(0), hit1_(0), hit2_(0), hit3_(0),
	  hitBoth23_(0), hitOnlyEj_(0), hitOnly1_(0), hitOnly2_(0), hitOnly3_(0),
	  hitOnly12_(0), hitOnly23_(0), hitOnly13_(0), hitOnly123_(0),
	  onePartHits_(0), twoPartHits_(0), threePartHits_(0), fourPartHits_(0),
	  detectorHits_(SABRE_Array.size()),
	  beamspot_(beamspot)
	{

	}

void det4mc::Run(std::ifstream& infile, std::ofstream& outfile){

	// std::cerr << "SABRE_Array_ size: " << SABRE_Array_.size()
	// 		  << ", SABREARRAY_EnergyResolutionModels_ size: " << SABREARRAY_EnergyResolutionModels_.size()
	// 		  << std::endl;


	double e1, theta1, phi1, e2, theta2, phi2, e3, theta3, phi3, e4, theta4, phi4;
	while(infile >> e1 >> theta1 >> phi1 >> e2 >> theta2 >> phi2 >> e3 >> theta3 >> phi3 >> e4 >> theta4 >> phi4){



		nevents_ += 1;
		bool detectedEj = false;	//kin4mc particle 1 == ejectile
		bool detected1 = false;		//kin4mc particle 2 == breakup 1 (first decay particle)
		bool detected2 = false;		//kin4mc particle 3 == breakup 2 (second decay particle)
		bool detected3 = false;		//kin4mc particle 4 == breakup 3 (third decay particle -- left over)

		//get reaction origin based on beamspot:
		Vec3 reactionOrigin = beamspot_->GeneratePoint();
		//Vec3 reactionOrigin = {0.,0.,0.};

		outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\t" << e3 << "\t" << theta3 << "\t" << phi3 << "\t" << e4 << "\t" << theta4 << "\t" << phi4 << std::endl;

		if(nevents_%50000==0) ConsoleColorizer::PrintBlue("Processed " + std::to_string(nevents_) + " events...\n");


		for(size_t i=0; i<SABRE_Array_.size(); i++){
			/*///////////////////////////////////////////////
			//  check if particle 1 (ejectile) is detected //
			///////////////////////////////////////////////*/
			double smearedERing = 0., smearedEWedge = 0.;



			std::pair<int,int> hit1_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD,reactionOrigin);

			if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detectedEj){
				//std::cout << "here1" << std::endl;
				//apply target energy loss to e1:
				double e1_aftertarget = targetLoss_par1_->ApplyEnergyLoss(e1, theta1);

				//apply dead layer energy loss to e1_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta1*DEG2RAD,phi1*DEG2RAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e1_afterDeadLayer = deadLayerLoss_par1_->ApplyEnergyLoss(e1_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit1_rw.first, e1_afterDeadLayer, smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit1_rw.second, e1_afterDeadLayer, smearedEWedge)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
					outfile << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detectedEj = true;
					hitej_ += 1;
					detectorHits_[i] += 1;
				}
			}


			/*//////////////////////////////////////////
			//  check if particle 2 (bu1) is detected //
			//////////////////////////////////////////*/
			smearedERing = 0.;
			smearedEWedge = 0.;

			std::pair<int,int> hit2_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta2*DEG2RAD,phi2*DEG2RAD,reactionOrigin);

			if(hit2_rw.first != -1 && hit2_rw.second != -1 && !detected1){
				//std::cout << "here2" << std::endl;
				//apply energy loss to e1:
				double e2_aftertarget = targetLoss_par2_->ApplyEnergyLoss(e2, theta2);

				//apply dead layer energy loss to e2_aftertarget
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta2*DEG2RAD,phi2*DEG2RAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal*(1/normal.Mag());
				double e2_afterDeadLayer = deadLayerLoss_par2_->ApplyEnergyLoss(e2_aftertarget, trajectory, normal);

				//std::cout << "here2 again" << std::endl;

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit2_rw.first,e2_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit2_rw.second,e2_afterDeadLayer,smearedEWedge)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit2_rw.first,hit2_rw.second);
					outfile << 200+i << "\t" << hit2_rw.first << "\t" << hit2_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected1 = true;
					hit1_ += 1;
					detectorHits_[i] += 1;
					//std::cout << "detectorHits_ incremented in here2" << std::endl;
				}

				//std::cout << "here2 again again" << std::endl;
			}

			/*//////////////////////////////////////////
			//  check if particle 3 (bu2) is detected //
			//////////////////////////////////////////*/
			smearedERing = 0.;
			smearedEWedge = 0.;

			std::pair<int,int> hit3_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta3*DEG2RAD,phi3*DEG2RAD,reactionOrigin);

			if(hit3_rw.first != -1 && hit3_rw.second != -1 && !detected2){
				//std::cout << "here3" << std::endl;
				//apply energy loss to e3:
				double e3_aftertarget = targetLoss_par3_->ApplyEnergyLoss(e3, theta3);

				//apply dead layer energy loss to e3_aftertarget:
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta3*DEG2RAD,phi3*DEG2RAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal * (1/normal.Mag());
				double e3_afterDeadLayer = deadLayerLoss_par3_->ApplyEnergyLoss(e3_aftertarget, trajectory, normal);

				//std::cout << "here3 again" << std::endl;

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit3_rw.first,e3_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit3_rw.second,e3_afterDeadLayer,smearedEWedge)){
					//std::cout << "here3 again again 1" << std::endl;
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit3_rw.first, hit3_rw.second);
					//std::cout << "here3 again again 2" << std::endl;
					outfile << 300+i << "\t" << hit3_rw.first << "\t" << hit3_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					//std::cout << "here3 again again 3" << std::endl;
					detected2 = true;
					//std::cout << "here3 again again 4" << std::endl;
					hit2_ += 1;
					//std::cout << "here3 again again 5" << std::endl;
					detectorHits_[i] += 1;
					//std::cout << "here3 again again 6" << std::endl;
				}

				//std::cout << "here3 again again again" << std::endl;

			}

			/*//////////////////////////////////////////
			//  check if particle 4 (bu2) is detected //
			//////////////////////////////////////////*/
			smearedERing = 0.;
			smearedEWedge = 0.;

			std::pair<int,int> hit4_rw = SABRE_Array_[i]->GetOffsetTrajectoryRingWedge(theta4*DEG2RAD,phi4*DEG2RAD,reactionOrigin);

			if(hit4_rw.first != -1 && hit4_rw.second != -1 && !detected3){
				//std::cout << "here4" << std::endl;
				//apply energy loss to e4:
				double e4_aftertarget = targetLoss_par4_->ApplyEnergyLoss(e4, theta4);

				//apply dead layer energy loss to e3_aftertarget
				Vec3 trajectory;
				trajectory.SetVectorSpherical(1,theta4*DEG2RAD,phi4*DEG2RAD);
				Vec3 normal = SABRE_Array_[i]->GetNormTilted();
				normal = normal * (1/normal.Mag());
				double e4_afterDeadLayer = deadLayerLoss_par4_->ApplyEnergyLoss(e4_aftertarget, trajectory, normal);

				if(SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInRing(hit4_rw.first,e4_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels_[i]->detectEnergyInWedge(hit4_rw.second,e4_afterDeadLayer,smearedEWedge)){
					Vec3 localCoords = SABRE_Array_[i]->GetHitCoordinatesRandomWiggle(hit4_rw.first,hit4_rw.second);
					outfile << 400+i << "\t" << hit4_rw.first << "\t" << hit4_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << std::endl;
					detected3 = true;
					hit3_ += 1;
					detectorHits_[i] += 1;
				}
			}

		//std::cout << "here end 1" << std::endl;

		}

		//std::cout << "here end 2" << std::endl;

		outfile << eoev << std::endl;

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
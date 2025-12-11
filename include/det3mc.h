#ifndef DET3MC_H
#define DET3MC_H

#include <string>
#include <vector>
#include <fstream>
#include "SABRE_Detector.h"
#include "SABRE_EnergyResolutionModel.h"
#include "TargetEnergyLoss.h"
#include "SABRE_DeadLayerModel.h"
#include "Vec3.h"
#include "TH1.h"
#include <cmath>
#include "Beamspot.h"
#include "RootWriter.h"
#include "plot3mc.h"
#include "TargetAngularStraggler.h"

class det3mc{
public:
	static constexpr double DEG2RAD = M_PI/180.;
	static constexpr double RAD2DEG = 180./M_PI;

	det3mc(std::vector<SABRE_Detector*>& SABRE_Array,
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
		   TargetAngularStraggler* straggler_par4);

	void Run(std::ifstream& infile, std::ofstream& outfile, RootWriter* RootWriter, plot3mc* RootPlotter, bool targetStraggle1, bool targetStraggle2, bool targetStraggle3, bool targetStraggle4);

	//after Run(), these functions can be called to query for statistics:
	long GetNumEvents() const;
	long GetHit1() const;
	long GetHit3() const;
	long GetHit4() const;
	long GetHitBoth34() const;
	long GetHitOnly3() const;
	long GetHitOnly4() const;
	long GetOnePartHits() const;
	long GetTwoPartHits() const;
	long GetThreePartHits() const;
	const std::vector<long>& GetDetectorHits() const;

	std::string GetToString_TargetLoss1() { return targetLoss_par1_->ToString(); }
	std::string GetToString_TargetLoss2() { return targetLoss_par2_->ToString(); }
	std::string GetToString_TargetLoss3() { return targetLoss_par3_->ToString(); }
	std::string GetToString_TargetLoss4() { return targetLoss_par4_->ToString(); }

	std::string GetToString_DeadLayerLoss1() { return deadLayerLoss_par1_->ToString(); }
	std::string GetToString_DeadLayerLoss2() { return deadLayerLoss_par2_->ToString(); }
	std::string GetToString_DeadLayerLoss3() { return deadLayerLoss_par3_->ToString(); }
	std::string GetToString_DeadLayerLoss4() { return deadLayerLoss_par4_->ToString(); }

private:
	std::vector<SABRE_Detector*> SABRE_Array_;
	std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels_;
	// TargetEnergyLoss* targetLoss_;
	// SABRE_DeadLayerModel* deadLayerLoss_;

	//	For 3-body final state:		A(a,b)B, B->c+d
	TargetEnergyLoss* targetLoss_par1_;//ejectile from kin3mc (b)
	TargetEnergyLoss* targetLoss_par2_;//recoil from kin3mc (B) - not detected!
	TargetEnergyLoss* targetLoss_par3_;//breakup from kin3mc (c)
	TargetEnergyLoss* targetLoss_par4_;//final daughter from kin3mc (d)

	SABRE_DeadLayerModel* deadLayerLoss_par1_;//ejectile from kin3mc (b)
	SABRE_DeadLayerModel* deadLayerLoss_par2_;//recoil from kin3mc (B) - not detected!
	SABRE_DeadLayerModel* deadLayerLoss_par3_;//breakup from kin3mc (c)
	SABRE_DeadLayerModel* deadLayerLoss_par4_;//final daughter from kin3mc (d)

	TargetAngularStraggler* straggler_par1_;
	TargetAngularStraggler* straggler_par2_;
	TargetAngularStraggler* straggler_par3_;
	TargetAngularStraggler* straggler_par4_;

	//counters for statistics:
	long nevents_;
	long hit1_;
	long hit3_;
	long hit4_;
	long hitBoth34_;
	long hitOnly3_;
	long hitOnly4_;
	long onePartHits_;
	long twoPartHits_;
	long threePartHits_;

	std::vector<long> detectorHits_;

	Beamspot* beamspot_;

	static const std::pair<int, int> offsets[];
};

#endif //DET3MC_H
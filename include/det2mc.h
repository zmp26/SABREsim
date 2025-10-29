#ifndef DET2MC_H
#define DET2MC_H

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

class det2mc {
public:
	static constexpr double DEG2RAD = M_PI/180.;
	static constexpr double RAD2DEG = 180./M_PI;

	det2mc(std::vector<SABRE_Detector*>& SABRE_Array,
		   std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels,
		   TargetEnergyLoss* targetLoss_par1,
		   TargetEnergyLoss* targetLoss_par2,
		   SABRE_DeadLayerModel* deadLayerLoss_par1,
		   SABRE_DeadLayerModel* deadLayerLoss_par2,
		   Beamspot* beamspot);

	void Run(std::ifstream& infile, std::ofstream& outfile);

	//after Run(), these functions may be called to query for statistics:
	long GetNumEvents() const;
	long GetHit1() const;
	long GetHit2() const;
	long GetHitBoth() const;
	long GetHit1Only() const;
	long GetHit2Only() const;
	const std::vector<long>& GetDetectorHits() const;

	std::string GetToString_TargetLoss1() { return targetLoss_par1_->ToString(); }
	std::string GetToString_TargetLoss2() { return targetLoss_par2_->ToString(); }

	std::string GetToString_DeadLayerLoss1() { return deadLayerLoss_par1_->ToString(); }
	std::string GetToString_DeadLayerLoss2() { return deadLayerLoss_par2_->ToString(); }

private:
	std::vector<SABRE_Detector*>& SABRE_Array_;
	std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels_;
	//TargetEnergyLoss* targetLoss_;
	//SABRE_DeadLayerModel* deadLayerLoss_;

	//	For 2-body reaction:		A(a,b)B
	TargetEnergyLoss* targetLoss_par1_;//ejectile from kin2mc (b)
	TargetEnergyLoss* targetLoss_par2_;//recoil from kin2mc	  (B)

	SABRE_DeadLayerModel* deadLayerLoss_par1_;//ejectile from kin2mc (b)
	SABRE_DeadLayerModel* deadLayerLoss_par2_;//recoil from kin2mc	  (B)

	//counters for statistics
	long nevents_;
	long hit1_;
	long hit2_;
	long hitBoth_;
	long hit1Only_;
	long hit2Only_;
	std::vector<long> detectorHits_;

	Beamspot* beamspot_;
};

#endif //DET2MC_H
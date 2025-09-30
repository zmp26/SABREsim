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

class det3mc{
public:
	static constexpr double DEG2RAD = M_PI/180.;
	static constexpr double RAD2DEG = 180./M_PI;

	det3mc(std::vector<SABRE_Detector*>& SABRE_Array,
		   std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels,
		   TargetEnergyLoss* targetLoss,
		   SABRE_DeadLayerModel* deadLayerLoss);

	void Run(std::ifstream& infile, std::ofstream& outfile);

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

private:
	std::vector<SABRE_Detector*> SABRE_Array_;
	std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels_;
	TargetEnergyLoss* targetLoss_;
	SABRE_DeadLayerModel* deadLayerLoss_;

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
};

#endif //DET3MC_H
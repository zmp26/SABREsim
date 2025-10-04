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
		   TargetEnergyLoss* targetLoss,
		   SABRE_DeadLayerModel* deadLayerLoss,
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

private:
	std::vector<SABRE_Detector*>& SABRE_Array_;
	std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels_;
	TargetEnergyLoss* targetLoss_;
	SABRE_DeadLayerModel* deadLayerLoss_;

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
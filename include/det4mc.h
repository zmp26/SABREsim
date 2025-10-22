#ifndef DET4MC_H
#define DET4MC_H

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

class det4mc{
public:
	static constexpr double DEG2RAD = M_PI/180.;
	static constexpr double RAD2DEG = 180./M_PI;

	det4mc(std::vector<SABRE_Detector*>& SABRE_Array,
		   std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels,
		   TargetEnergyLoss* targetLossEj, TargetEnergyLoss* targetLoss1, TargetEnergyLoss* targetLoss2, TargetEnergyLoss* targetLoss3,
		   SABRE_DeadLayerModel* deadLayerLossEj, SABRE_DeadLayerModel* deadLayerLoss1, SABRE_DeadLayerModel* deadLayerLoss2, SABRE_DeadLayerModel* deadLayerLoss3,
		   Beamspot* beamspot);

	void Run(std::ifstream& infile, std::ofstream& outfile);

	/*
		A(a,b)B
	      	   \
	            B → D + c
	                 \
	                  D → d + e
	*/

	//after Run(), these functions can be called to query for statistics:
	long GetNumEvents() const;
	long GetHitEj() const;									//ejectile 						(b)
	long GetHit1() const;									//at least decay 1 				(c)
	long GetHit2() const;									//at least decay 2 				(d)
	long GetHit3() const;									//at least decay 3 				(e)
	long GetHitBoth23() const;								//at least decay2, decay 3 		(d,e)
	long GetHitOnlyEj() const;								//ejectile only 				(b)
	long GetHitOnly1() const;								//decay 1 only 					(c)
	long GetHitOnly2() const;								//decay 2 only 					(d)
	long GetHitOnly3() const;								//decay 3 only 					(e)
	long GetHitOnly12() const;								//decay 1&2 only 				(c,d)
	long GetHitOnly23() const;								//decay 2&3 only 				(d,e)
	long GetHitOnly13() const;								//decay 1&3 only 				(c,e)
	long GetHitOnly123() const;								//decay 1,2,3 only 				(c,d,e)
	long GetOnePartHits() const;							//total events w/ 1 particle
	long GetTwoPartHits() const;							//total events w/ 2 particles
	long GetThreePartHits() const;							//total events w/ 3 particles
	long GetFourPartHits() const;							//total events w/ 4 particles
	const std::vector<long>& GetDetectorHits() const;		//SABRE array total hits/det

private:
	std::vector<SABRE_Detector*> SABRE_Array_;
	std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels_;
	TargetEnergyLoss* targetLossEj_;
	TargetEnergyLoss* targetLoss1_;
	TargetEnergyLoss* targetLoss2_;
	TargetEnergyLoss* targetLoss3_;
	SABRE_DeadLayerModel* deadLayerLossEj_;
	SABRE_DeadLayerModel* deadLayerLoss1_;
	SABRE_DeadLayerModel* deadLayerLoss2_;
	SABRE_DeadLayerModel* deadLayerLoss3_;

	//counters for statistics:
	long nevents_;
	long hitej_;
	long hit1_;
	long hit2_;
	long hit3_;
	long hitBoth23_;
	long hitOnlyEj_;
	long hitOnly1_;
	long hitOnly2_;
	long hitOnly3_;
	long hitOnly12_;//new
	long hitOnly23_;
	long hitOnly13_;//new
	long hitOnly123_;
	long onePartHits_;
	long twoPartHits_;
	long threePartHits_;
	long fourPartHits_;

	std::vector<long> detectorHits_;

	Beamspot* beamspot_;

};

#endif//DET4MC_H
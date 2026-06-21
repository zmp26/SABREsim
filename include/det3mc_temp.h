#ifndef DET3MC_TEMP_H
#define DET3MC_TEMP_H

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "structs.h"
#include "SABRE_Array.h"
#include "plot3mc.h"
#include "EventRecorder.h"
#include "SimConfig.h"
#include "Beamspot.h"
#include "SPS_Aperture.h"
#include "SABRE_EnergyResolutionModel.h"

class det3mc_temp {
public:
	det3mc_temp(SABRE_Array* SABRE_Array_, SPS_Aperture* SPS_Aperture_, 
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


	void Run(std::ifstream& infile, std::ofstream& outfile, EventRecorder* EventRecorder, plot3mc* RootPlotter, SimConfig* config);

	/*
		A(a,b)B
	      	   \
	            B → c + d
	*/

	//after Run(), these functions can be called to query for statistics:
	long GetNumEvents() const { return nevents_; };								//number of events considered by simulation (e.g. does not count events w/o EjInSPS if coincidence = true)
	long GetNumKinEvents() const { return nkinevents_; };							//number of kinematic events in input file
	long GetNoHitEjSPS() const { return nohitejSPS_; };								//number of input events where ejectile does not intersect SPS region
	long GetHitEjSPS() const { return hitejSPS_; };								//ejectile in SPS
	long GetHit1() const { return hit1_; };									//ejectile in SABRE				(b)
	long GetHit2() const { return hit2_; };									//at least decay 1 				(c)
	long GetHit3() const { return hit3_; };									//at least decay 2				(d)
	// long GetHit3() const { return hit3_; };									//at least decay 3 				(e)
	long GetHitBoth23() const { return hitBoth23_; };								//at least decay2, decay 3 		(d,e)
	//long GetHitOnlyEj() const { return hitOnlyEj_; };								//ejectile only 				(b)
	long GetHitOnly1() const { return hitOnly1_; };								//decay 1 only 					(c)
	long GetHitOnly2() const { return hitOnly2_; };								//decay 2 only 					(d)
	long GetHitOnly3() const { return hitOnly3_; };								//decay 3 only 					(e)
	long GetHitOnly12() const { return hitOnly12_; };								//decay 1&2 only 				(c,d)
	long GetHitOnly23() const { return hitOnly23_; };								//decay 2&3 only 				(d,e)
	long GetHitOnly13() const { return hitOnly13_; };								//decay 1&3 only 				(c,e)
	long GetHitOnly123() const { return hitOnly123_; };								//decay 1,2,3 only 				(c,d,e)
	long GetOnePartHits() const { return onePartHits_; };							//total events w/ 1 particle
	long GetTwoPartHits() const { return twoPartHits_; };							//total events w/ 2 particles
	long GetThreePartHits() const { return threePartHits_; };							//total events w/ 3 particles
	const std::vector<long>& GetDetectorHits() const { return detectorHits_;};

	std::string GetToString_TargetLoss1() { return targetLoss_par1_->ToString(); }
	std::string GetToString_TargetLoss2() { return targetLoss_par2_->ToString(); }
	std::string GetToString_TargetLoss3() { return targetLoss_par3_->ToString(); }
	std::string GetToString_TargetLoss4() { return targetLoss_par4_->ToString(); }

	std::string GetToString_DeadLayerLoss1() { return deadLayerLoss_par1_->ToString(); }
	std::string GetToString_DeadLayerLoss2() { return deadLayerLoss_par2_->ToString(); }
	std::string GetToString_DeadLayerLoss3() { return deadLayerLoss_par3_->ToString(); }
	std::string GetToString_DeadLayerLoss4() { return deadLayerLoss_par4_->ToString(); }


private:
	//long GetHitOnlyEj() const { return hitOnlyEj_; };								//ejectile only 				(b)

	double DEGRAD = M_PI/180.;
	double RADDEG = 180./M_PI;

	SABRE_Array* SABRE_Array_;
	SPS_Aperture* SPS_Aperture_;
	std::vector<SABRE_EnergyResolutionModel*>& SABREARRAY_EnergyResolutionModels_;

	//For 4-body final state3		A(a,b)B, B->D+c, D->d+e   (end w/ b, c, d, e in final state)
	TargetEnergyLoss* targetLoss_par1_;	//ejectile from kin4mc (b)
	TargetEnergyLoss* targetLoss_par2_; //breakup1 from kin4mc (c)
	TargetEnergyLoss* targetLoss_par3_; //breakup2 from kin4mc (d)
	TargetEnergyLoss* targetLoss_par4_; //breakup3/final daughter from kin4mc (e)

	SABRE_DeadLayerModel* deadLayerLoss_par1_; //ejectile from kin4mc (b)
	SABRE_DeadLayerModel* deadLayerLoss_par2_; //breakup1 from kin4mc (c)
	SABRE_DeadLayerModel* deadLayerLoss_par3_; //breakup2 from kin4mc (d)
	SABRE_DeadLayerModel* deadLayerLoss_par4_; //breakup3/final daughter from kin4mc (e)

	TargetAngularStraggler* straggler_par1_;
	TargetAngularStraggler* straggler_par2_;
	TargetAngularStraggler* straggler_par3_;
	TargetAngularStraggler* straggler_par4_;

	//counters for statistics:
	long nkinevents_;
	long nevents_;
	long hitejSPS_;
	long nohitejSPS_;
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

	std::vector<long> detectorHits_;

	Beamspot* beamspot_;

	static const std::pair<int, int> offsets[];

	bool ProcessParticle(Particle& p, const Vec3& origin, EventRecorder* rec, plot3mc* plotter, SimConfig* config, std::ostringstream& ss);

};

#endif//DET3MC_TEMP_H
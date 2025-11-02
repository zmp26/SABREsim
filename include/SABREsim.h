#ifndef SABRESIM_H
#define SABRESIM_H

#include <string>
#include <vector>
#include <fstream>
#include "TH1.h"
#include "SABRE_Detector.h"
#include "SABRE_EnergyResolutionModel.h"
#include "TargetEnergyLoss.h"
#include "SABRE_DeadLayerModel.h"
#include "det2mc.h"
#include "det3mc.h"
#include "det4mc.h"
#include "RootWriter.h"

class SABREsim {
public:
	SABREsim(int kinX,
			 const std::string& kinInputFilename,
			 const std::string& detOutputFilename);

	~SABREsim();

	void Run();

private:
	static constexpr double DEG2RAD = M_PI / 180.0;
	static constexpr double RAD2DEG = 180.0 / M_PI;

	int kinX_;
	std::string kinInputFilename_;
	std::string detOutputFilename_;

	std::vector<SABRE_Detector*> SABRE_Array_;
	std::vector<SABRE_EnergyResolutionModel*> SABREARRAY_EnergyResolutionModels_;

	//---------------------------TARGET ENERGY LOSSES------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	/*
			BELOW:
			These hold the information for all possible energy losses we have made. Easy to just read them in since they are ultimately just TF1 plus a few extra pieces of info

			These are what the targetLoss_par1/2/3/4_ should be set to before being passed to det2/3/4mc classes

	*/

	TargetEnergyLoss* targetLoss_6Li_in_LiF_;
	TargetEnergyLoss* targetLoss_alpha_in_LiF_;
	TargetEnergyLoss* targetLoss_deuteron_in_LiF_;
	TargetEnergyLoss* targetLoss_none_;

	TargetEnergyLoss* targetLoss_alpha_in_10B_;//NEED TO MAKE THIS CONFIG FILE AND FIT STILL!!!!!!!!!!!!!!!*****************************************************************<=======
	TargetEnergyLoss* targetLoss_proton_in_10B_;//NEED TO MAKE THIS CONFIG FILE AND FIT STILL!!!!!!!!!!!!!!!*****************************************************************<=======

	//TargetEnergyLoss* targetLoss_; //pointer to current TargetEnergyLoss in use!

	/*
			BELOW:
			Current Target Energy Loss pointers -- these are the target energy loss actually passed to det2/3/4mc

			For det2mc:
				targetLoss_par1_ = ejectile energy loss -> this can be set to targetLoss_none_ as we are ignoring the ejectile if using SPS acceptances. However, the handling is general and this *can* also be something else for singles data, etc.
				targetLoss_par2_ = recoil energy loss -> this can be set to whichever target loss for recoil from kin2mc file
				targetLoss_par3_ = targetLoss_none_ as there is no particle 3
				targetLoss_par4_ = targetLoss_none_ as there is no particle 4

			For det3mc:
				targetLoss_par1_ = ejectile energy loss -> this can be set to targetLoss_none_ as we are ignoring the ejectile if using SPS accepatances. However, the handling is general and this *can* also be something else for singles data, etc.
				targetLoss_par2_ = recoil energy loss -> this can be set to targetLoss_none_ as we are not detecting the recoil, only its break ups
				targetLoss_par3_ = breakup1 energy loss -> this can be set to targetLoss for whatever the specified breakup particle is in kin3mc
				targetLoss_par4_ = breakup2/daughter energy loss -> this can be set to targetLoss for whatever the daughter particle is in kin3mc

			for det4mc:
				targetLoss_par1_ = ejectile energy loss -> this can be set to targetLoss_none_ as we are ignoring the ejectile if using SPS accepatances. However, the handling is general and this *can* also be something else for singles data, etc.
				targetLoss_par2_ = breakup1 energy loss -> this can be set to appropriate energy loss for first designated breakup particle from kin4mc
				targetLoss_par3_ = breakup2 energy loss -> this can be set to appropriate energy loss for second designated breakup particle from kin4mc
				targetLoss_par4_ = breakup3/final daughter energy loss -> this can be set to appropriate energy loss for final daughter from kin4mc
	*/

	TargetEnergyLoss* targetLoss_par1_;
	TargetEnergyLoss* targetLoss_par2_;
	TargetEnergyLoss* targetLoss_par3_;
	TargetEnergyLoss* targetLoss_par4_;


	//---------------------------DEAD LAYER ENERGY LOSSES-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	/*
			BELOW:
			These hold the information for all possible dead layer energy losses we have made. Easy to just read them in since they are ultimately just TF1 plus a few extra pieces of info

			These are what the deadLayerLoss_par1/2/3/4_ should be set to before being passed to det2/3/4mc classes

	*/

	SABRE_DeadLayerModel* deadLayerLoss_6Li_;
	SABRE_DeadLayerModel* deadLayerLoss_alpha_;
	SABRE_DeadLayerModel* deadLayerLoss_deuteron_;
	SABRE_DeadLayerModel* deadLayerLoss_proton_;//NEED TO MAKE THIS CONFIG FILE AND FIT STILL!!!!!!!!!!!!!!!*****************************************************************<=======
	SABRE_DeadLayerModel* deadLayerLoss_none_;

	//SABRE_DeadLayerModel* deadLayerLoss_; //pointer to current SABRE_DeadLayerModel in use!

	/*
			BELOW:
			Current Dead Layer Loss pointers -- these are the dead layer energy loss actually passed to det2/3/4mc

			for det2mc:
				deadLayerLoss_par1_ = ejectile dead layer loss - set to none unless singles data, etc.
				deadLayerLoss_par2_ = recoil dead layer loss - set to appropriate value for recoil from kin2mc
				deadLayerLoss_par3_ = none, not used
				deadLayerLoss_par4_ = none, not used

			for det3mc:
				deadLayerLoss_par1_ = ejectile dead layer loss - set to none unless singles data, etc.
				deadLayerLoss_par2_ = recoil dead layer loss - set to none since we do not measure recoil
				deadLayerLoss_par3_ = breakup1 dead layer loss  - set to appropriate value for breakup from kin3mc
				deadLayerLoss_par4_ = breakup2/deaughter dead layer loss - set to appropriate value for breakup2/daughter from kin3mc

			for det4mc:
				deadLayerLoss_par1_ = ejectile dead layer loss - set to none unless singles data, etc.
				deadLayerLoss_par2_ = breakup1 dead layer loss - set to appropriate value for breakup1 from kin4mc
				deadLayerLoss_par3_ = breakup2 dead layer loss - set to appropraite value for breakup2 from kin4mc
				deadLayerLoss_par4_ = brekaup3/final daughter dead layer loss - set to appropriate value for breakup3/final daughter from kin4mc
	*/

	SABRE_DeadLayerModel* deadLayerLoss_par1_;
	SABRE_DeadLayerModel* deadLayerLoss_par2_;
	SABRE_DeadLayerModel* deadLayerLoss_par3_;
	SABRE_DeadLayerModel* deadLayerLoss_par4_;

	BeamProfile* profile_;
	Beamspot* beamspot_;

	RootWriter* RootWriter_;

	long nevents_;
	std::vector<long> detectorHits_;
	long hit1_, hit2_, hit3_, hit4_, hit23_, hitej_;
	long hit34_, hitOnly1_, hitOnly2_, hitOnly3_, hitOnly4_, hitOnlyEj_, hitOnly12_, hitOnly23_, hitOnly13_, hitOnly123_;
	long hit1Only_, hit2Only_, hitBoth_;
	long onePartHits_, twoPartHits_, threePartHits_, fourPartHits_;

	void InitializeDetectors(bool WriteCornersToFile=false);
	bool InitializeModels();
	void InitializeBeamspot();
	void CleanUp();

	void Simulate2body(std::ifstream& infile, std::ofstream& outfile);
	void Simulate3body(std::ifstream& infile, std::ofstream& outfile);
	void Simulate4body(std::ifstream& infile, std::ofstream& outfile);

	void PrintSummary() const;
};

#endif //SABRESIM_H
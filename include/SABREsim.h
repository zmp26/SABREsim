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

	TargetEnergyLoss* targetLoss_6Li_in_LiF_;
	TargetEnergyLoss* targetLoss_alpha_in_LiF_;
	TargetEnergyLoss* targetLoss_deuteron_in_LiF_;

	TargetEnergyLoss* targetLoss_; //pointer to current TargetEnergyLoss in use!

	SABRE_DeadLayerModel* deadLayerLoss_6Li_;
	SABRE_DeadLayerModel* deadLayerLoss_alpha_;
	SABRE_DeadLayerModel* deadLayerLoss_deuteron_;

	SABRE_DeadLayerModel* deadLayerLoss_; //pointer to current SABRE_DeadLayerModel in use!

	BeamProfile* profile_;
	Beamspot* beamspot_;

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
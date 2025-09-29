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

	long nevents_;
	std::vector<long> detectorHits_;
	long hit1_, hit2_, hit3_, hit4_;
	long hit34_, hitOnly3_, hitOnly4_;
	long hit1Only_, hit2Only_, hitBoth_;
	long onePartHits_, twoPartHits_, threePartHits_, fourPartHits_;

	void InitializeDetectors();
	bool InitializeModels();
	void CleanUp();

	void Simulate2body(std::ifstream& infile, std::ofstream& outfile);
	void Simulate3body(std::ifstream& infile, std::ofstream& outfile);
	void Simulate4body(std::ifstream& infile, std::ofstream& outfile);

	void PrintSummary() const;
};

#endif //SABRESIM_H
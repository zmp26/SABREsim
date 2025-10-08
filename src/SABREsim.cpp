#include "SABREsim.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include "ConsoleColorizer.h"
#include "Beamspot.h"
#include "UniformProfile.h"
#include "GaussianProfile.h"

static const std::pair<int,int> offsets[] = {
	{112,40}, {96,32}, {80,16}, {64,24}, {48,0}
};

SABREsim::SABREsim(int kinX,
				   const std::string& kinInputFilename,
				   const std::string& detOutputFilename)
	: kinX_(kinX),
	  kinInputFilename_(kinInputFilename),
	  detOutputFilename_(detOutputFilename),
	  targetLoss_6Li_in_LiF_(nullptr),
	  targetLoss_alpha_in_LiF_(nullptr),
	  targetLoss_deuteron_in_LiF_(nullptr),
	  deadLayerLoss_6Li_(nullptr),
	  deadLayerLoss_alpha_(nullptr),
	  deadLayerLoss_deuteron_(nullptr),
	  nevents_(0),
	  detectorHits_(5,0),
	  hit1_(0), hit2_(0), hit3_(0), hit4_(0),
      hit34_(0), hitOnly3_(0), hitOnly4_(0),
      hit1Only_(0), hit2Only_(0), hitBoth_(0),
      onePartHits_(0), twoPartHits_(0), threePartHits_(0), fourPartHits_(0)
{
	srand(static_cast<unsigned>(time(nullptr)));
}

SABREsim::~SABREsim(){
	CleanUp();
}

void SABREsim::CleanUp(){
	for (auto d : SABRE_Array_) delete d;
	SABRE_Array_.clear();

	for(auto m : SABREARRAY_EnergyResolutionModels_) delete m;
	SABREARRAY_EnergyResolutionModels_.clear();

	delete targetLoss_6Li_in_LiF_;
	delete targetLoss_alpha_in_LiF_;
	delete targetLoss_deuteron_in_LiF_;

	delete deadLayerLoss_6Li_;
	delete deadLayerLoss_alpha_;
	delete deadLayerLoss_deuteron_;

}

void SABREsim::InitializeDetectors(bool WriteCornersToFile){
	double INNER_R = 0.0326;
	double OUTER_R = 0.1351;
	double TILT = 40.0;
	double ZDIST = -0.1245;
	double PHI_COVERAGE = 54.4;
	std::vector<double> PHI = {306., 18., 234., 162., 90.};

	for(size_t i=0; i<PHI.size(); i++){
		SABRE_Detector* det = new SABRE_Detector(INNER_R, OUTER_R, PHI_COVERAGE*DEG2RAD, PHI[i]*DEG2RAD, TILT*DEG2RAD, ZDIST);
		SABRE_Array_.push_back(det);
		// std::cout << "Created SABRE_Detector[" << i << "] at phi = " << PHI[i]
		//           << " norm = (" << det->GetNormTilted.GetX()
		//           << ", " << det->GetNormTilted().GetY()
		//           << ", " << det->GetNormTilted().GetZ() << ")\n";
		Vec3 normtilted = det->GetNormTilted();
		TString outline = Form("Created SABRE_Detector[%zu] at phi = %f norm = (%f, %f, %f)\n",i,PHI[i],normtilted.GetX(),normtilted.GetY(),normtilted.GetZ());
		ConsoleColorizer::PrintGreen(outline.Data());

		if(WriteCornersToFile){

			TString cornerfilename = Form("config/corners/SABRE%zu_corners.txt",i);
			std::ofstream cornerfile(cornerfilename);

			det->WriteTransformedCorners(cornerfile);
		}

	}
	std::cout << "\n";

}

bool SABREsim::InitializeModels(){
	//energy resolution models
	for(size_t i=0; i<SABRE_Array_.size(); i++){
		SABRE_EnergyResolutionModel* m = new SABRE_EnergyResolutionModel(0.050, 0.100);
		SABREARRAY_EnergyResolutionModels_.push_back(m);
	}

	//load target energy losses
	targetLoss_6Li_in_LiF_ = TargetEnergyLoss::LoadFromConfigFile("config/TargetELoss_6Li_in_LiF.conf");
	if(!targetLoss_6Li_in_LiF_){
		ConsoleColorizer::PrintRed("Failed to load 6Li in LiF target loss config file\n");
		return false;
	}

	targetLoss_alpha_in_LiF_ = TargetEnergyLoss::LoadFromConfigFile("config/TargetELoss_alpha_in_LiF.conf");
	if(!targetLoss_alpha_in_LiF_){
		ConsoleColorizer::PrintRed("Failed to load alpha in LiF target loss config file\n");
		return false;
	}

	targetLoss_deuteron_in_LiF_ = TargetEnergyLoss::LoadFromConfigFile("config/TargetELoss_deuteron_in_LiF.conf");
	if(!targetLoss_deuteron_in_LiF_){
		ConsoleColorizer::PrintRed("Failed to load deuteron in LiF target loss config file\n");
		return false;
	}


	//load dead layer models
	deadLayerLoss_6Li_ = SABRE_DeadLayerModel::LoadFromConfigFile("config/DeadLayerELoss_6Li_in_Si.conf");
	if(!deadLayerLoss_6Li_){
		ConsoleColorizer::PrintRed("Failed to load 6Li dead layer config file\n");
		return false;
	}

	deadLayerLoss_alpha_ = SABRE_DeadLayerModel::LoadFromConfigFile("config/DeadLayerELoss_alpha_in_Si.conf");
	if(!deadLayerLoss_alpha_){
		ConsoleColorizer::PrintRed("Failed to load alpha dead layer config file\n");
		return false;
	}

	deadLayerLoss_deuteron_ = SABRE_DeadLayerModel::LoadFromConfigFile("config/DeadLayerELoss_deuteron_in_Si.conf");
	if(!deadLayerLoss_deuteron_){
		ConsoleColorizer::PrintRed("Failed to load deuteron dead layer config file\n");
		return false;
	}

	//SET CURRENT TARGETLOSS AND DEADLAYERLOSS HERE:
	//EVENTUALLY, THIS SHOULD BE DETERMINED BY A FILE READ IN TO AVOID HAVING TO REMAKE EVERY TIME
	targetLoss_ = targetLoss_6Li_in_LiF_;
	deadLayerLoss_ = deadLayerLoss_6Li_;

	return true;
}

void SABREsim::InitializeBeamspot(){
	//eventually interface with file read in so no need to remake when changing!

	profile_ = new GaussianProfile(0.002, 0.002);//meters

	beamspot_ = new Beamspot();
	beamspot_->SetProfile(profile_);
	beamspot_->SetBeamAxisOffset(0.,0.);//aligned along z axis (no lateral offset)
}

void SABREsim::Run(){
	std::ifstream infile(kinInputFilename_);
	if(!infile){
		ConsoleColorizer::PrintRed(("Error opening input file: " + kinInputFilename_).c_str());
		return;
	}

	std::ofstream outfile(detOutputFilename_);
	if(!outfile){
		ConsoleColorizer::PrintRed(("Error opening output file: " + detOutputFilename_).c_str());
		return;
	}

	InitializeDetectors(true);
	InitializeBeamspot();
	if(!InitializeModels()){
		ConsoleColorizer::PrintRed("Model initializiation failed. Exiting...\n");
		return;
	}

	if(kinX_ == 2){
		Simulate2body(infile, outfile);
	} else if(kinX_ == 3){
		Simulate3body(infile,outfile);
	} else if(kinX_ == 4){
		Simulate4body(infile,outfile);
	} else if(kinX_ == 0){ 
		//temporarily test Beamspot here!

		/* Establish beamspots here */
		Beamspot uniformbeam, gaussianbeam;

		uniformbeam.SetProfile(new UniformProfile(0.01,0.01));
		uniformbeam.SetBeamAxisOffset(0.,0.);

		gaussianbeam.SetProfile(new GaussianProfile(0.01, 0.01));
		gaussianbeam.SetBeamAxisOffset(0.,0.);

		TH2D* uniformbeam_xy = new TH2D("uniformbeam_xy","uniformbeam_xy", 1000, -0.05, 0.05, 1000, -0.05, 0.05);
		TH2D* gaussianbeam_xy = new TH2D("gaussianbeam_xy","gaussianbeam_xy", 1000, -0.05, 0.05, 1000, -0.05, 0.05);
		TH1D* uniformbeam_y = new TH1D("uniformbeam_y","uniformbeam_y",1000, -0.05, 0.05);
		TH1D* uniformbeam_x = new TH1D("uniformbeam_x","uniformbeam_x",1000, -0.05, 0.05);
		TH1D* gaussianbeam_x = new TH1D("gaussianbeam_x","gaussianbeam_x",1000, -0.05, 0.05);
		TH1D* gaussianbeam_y = new TH1D("gaussianbeam_y","gaussianbeam_y",1000, -0.05, 0.05);

		int numevents = 1000000; 
		for(int i=0; i<numevents; i++){
			if(i%50000 == 0) ConsoleColorizer::PrintBlue(std::string(Form("Processed %d events...\n",i)));

			/* Generate a point at z=0 */
			Vec3 uniformbeam_point = uniformbeam.GeneratePoint(0.);
			Vec3 gaussianbeam_point = gaussianbeam.GeneratePoint(0.);

			uniformbeam_xy->Fill(uniformbeam_point.GetX(),uniformbeam_point.GetY());
			gaussianbeam_xy->Fill(gaussianbeam_point.GetX(), gaussianbeam_point.GetY());
			uniformbeam_x->Fill(uniformbeam_point.GetX());
			uniformbeam_y->Fill(uniformbeam_point.GetY());
			gaussianbeam_x->Fill(gaussianbeam_point.GetX());
			gaussianbeam_y->Fill(gaussianbeam_point.GetY());
		}

		ConsoleColorizer::PrintBlue(std::string(Form("\n\nProcessed %d events!\n",numevents)));

		TCanvas* c1 = new TCanvas("c1","c1",800,600);
		uniformbeam_xy->Draw();
		c1->SaveAs("uniformbeam_xy.png");
		gaussianbeam_xy->Draw();
		c1->SaveAs("gaussianbeam_xy.png");

		//delete
		delete uniformbeam_xy;
		delete gaussianbeam_xy;
		delete uniformbeam_x;
		delete gaussianbeam_x;
		delete uniformbeam_y;
		delete gaussianbeam_y;
		delete c1;

	} else {
		//ConsoleColorizer::PrintRed(Form("Invalid kinX = %d\n",kinX_).Data());
		ConsoleColorizer::PrintRed("Invalid kinX = " + std::to_string(kinX_) + "\n");
	}

	infile.close();
	outfile.close();

	PrintSummary();
}

void SABREsim::Simulate2body(std::ifstream& infile, std::ofstream& outfile){
	//std::cout << "passing..." << std::endl;
	
	if(!infile.is_open()){
		ConsoleColorizer::PrintRed("Error: Cannot open input file!\n");
		return;
	}
	
	if(!outfile.is_open()){
		ConsoleColorizer::PrintRed("Error: Cannot open output file!\n");
		infile.close();
		return;
	}

	det2mc det2mcProcessor(SABRE_Array_, SABREARRAY_EnergyResolutionModels_, targetLoss_, deadLayerLoss_, beamspot_);

	det2mcProcessor.Run(infile,outfile);

	nevents_ = det2mcProcessor.GetNumEvents();
	detectorHits_ = det2mcProcessor.GetDetectorHits();
	hit1_ = det2mcProcessor.GetHit1();
	hit2_ = det2mcProcessor.GetHit2();
	hitBoth_ = det2mcProcessor.GetHitBoth();
	hit1Only_ = det2mcProcessor.GetHit1Only();
	hit2Only_ = det2mcProcessor.GetHit2Only();
}

void SABREsim::Simulate3body(std::ifstream& infile, std::ofstream& outfile){
	//std::cout << "passing..." << std::endl;

	if(!infile.is_open()){
		ConsoleColorizer::PrintRed("Error: Cannot open input file!\n");
		return;
	}

	if(!outfile.is_open()){
		ConsoleColorizer::PrintRed("Error: Cannot open output file!\n");
		return;
	}

	det3mc det3mcProcessor(SABRE_Array_, SABREARRAY_EnergyResolutionModels_, targetLoss_, deadLayerLoss_, beamspot_);

	det3mcProcessor.Run(infile, outfile);

	nevents_ = det3mcProcessor.GetNumEvents();
	detectorHits_ = det3mcProcessor.GetDetectorHits();
	hit1_ = det3mcProcessor.GetHit1();
	hit3_ = det3mcProcessor.GetHit3();
	hit4_ = det3mcProcessor.GetHit4();
	hit34_ = det3mcProcessor.GetHitBoth34();
	hitOnly3_ = det3mcProcessor.GetHitOnly3();
	hitOnly4_ = det3mcProcessor.GetHitOnly4();
	onePartHits_ = det3mcProcessor.GetOnePartHits();
	twoPartHits_ = det3mcProcessor.GetTwoPartHits();
	threePartHits_ = det3mcProcessor.GetThreePartHits();
}

void SABREsim::Simulate4body(std::ifstream& infile, std::ofstream& outfile){
	std::cout << "passing..." << std::endl;
}

void SABREsim::PrintSummary() const {
	std::cout << "\nProcessed " << nevents_ << " kin" << kinX_ << "mc events.\n";

	if (kinX_ == 2 || kinX_ == 3 || kinX_ == 4) {
		std::cout << "Events with ejectile in SABRE: " << hit1_
				  << " (" << float(hit1_)*100.0f / float(nevents_) << "%)\n";
	}
	if (kinX_ == 2) {
		std::cout << "Events with recoil in SABRE: " << hit2_
				  << " (" << float(hit2_)*100.0f / float(nevents_) << "%)\n";
		std::cout << "Only ejectile: " << hit1Only_
				  << " (" << float(hit1Only_)*100.0f / float(nevents_) << "%)\n";
		std::cout << "Only recoil: " << hit2Only_
				  << " (" << float(hit2Only_)*100.0f / float(nevents_) << "%)\n";
		std::cout << "Both: " << hitBoth_
				  << " (" << float(hitBoth_)*100.0f / float(nevents_) << "%)\n";
	}
	if (kinX_ == 3) {
		std::cout << "At least bu1: " << hit3_ << "\n";
		std::cout << "At least bu2: " << hit4_ << "\n";
		std::cout << "Only bu1: " << hitOnly3_
				  << " (" << float(hitOnly3_)*100.0f / float(nevents_) << "%)\n";
		std::cout << "Only bu2: " << hitOnly4_
				  << " (" << float(hitOnly4_)*100.0f / float(nevents_) << "%)\n";
		std::cout << "Both bu1 & bu2: " << hit34_
				  << " (" << float(hit34_)*100.0f / float(nevents_) << "%)\n";
	}
	if (kinX_ == 4) {
		std::cout << "1‑particle events: " << onePartHits_ << "\n";
		std::cout << "2‑particle events: " << twoPartHits_ << "\n";
		std::cout << "3‑particle events: " << threePartHits_ << "\n";
		std::cout << "4‑particle events: " << fourPartHits_ << "\n";
	}
	for (size_t i = 0; i < detectorHits_.size(); i++) {
		std::cout << "Detector_" << i << " total hits = " << detectorHits_[i] << "\n";
	}
	std::cout << "Output file: " << detOutputFilename_ << "\n\n";
}
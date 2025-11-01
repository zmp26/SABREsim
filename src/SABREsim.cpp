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
#include "FixedPointProfile.h"

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
	  targetLoss_none_(nullptr),
	  targetLoss_alpha_in_10B_(nullptr),
	  targetLoss_proton_in_10B_(nullptr),
	  deadLayerLoss_6Li_(nullptr),
	  deadLayerLoss_alpha_(nullptr),
	  deadLayerLoss_deuteron_(nullptr),
	  deadLayerLoss_proton_(nullptr),
	  deadLayerLoss_none_(nullptr),
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
			std::ofstream cornerfile_cartesian(cornerfilename);

			cornerfilename = Form("config/corners/SABRE%zu_corners_spherical.txt",i);
			std::ofstream cornerfile_spherical(cornerfilename);

			det->WriteTransformedCorners(cornerfile_cartesian);
			det->WriteTransformedCornersSpherical(cornerfile_spherical);

			cornerfile_cartesian.close();
			cornerfile_spherical.close();
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

	targetLoss_none_ = TargetEnergyLoss::LoadFromConfigFile("config/TargetELoss_none.conf");
	if(!targetLoss_none_){
		ConsoleColorizer::PrintRed("Failed to load none target loss config file\n");
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

	deadLayerLoss_none_ = SABRE_DeadLayerModel::LoadFromConfigFile("config/DeadLayerELoss_none.conf");
	if(!deadLayerLoss_none_){
		ConsoleColorizer::PrintRed("Failed to load none dead layer config\n");
		return false;
	}

	//SET CURRENT TARGETLOSS AND DEADLAYERLOSS HERE:
	//EVENTUALLY, THIS SHOULD BE DETERMINED BY A FILE READ IN TO AVOID HAVING TO REMAKE EVERY TIME
	// targetLoss_ = targetLoss_6Li_in_LiF_;
	// deadLayerLoss_ = deadLayerLoss_6Li_;

	return true;
}

void SABREsim::InitializeBeamspot(){
	//eventually interface with file read in so no need to remake when changing!

	//profile_ = new GaussianProfile(0.001, 0.001);//meters, 0.001m = 1mm
	//profile_ = new GaussianProfile(0.002, 0.002);//meters, 0.002m = 2mm
	//profile_ = new GaussianProfile(0.003, 0.003);//meters, 0.003m = 3mm
	//profile_ = new GaussianProfile(0.004, 0.004);
	//profile_ = new GaussianProfile(0.005, 0.005);
	//profile_ = new GaussianProfile(0.006, 0.006);
	//profile_ = new GaussianProfile(0.007, 0.007);
	//profile_ = new GaussianProfile(0.008,0.008);
	//profile_ = new GaussianProfile(0.009,0.009);
	//profile_ = new GaussianProfile(0.010,0.010);
	profile_ = new FixedPointProfile();

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

	InitializeDetectors();
	InitializeBeamspot();
	if(!InitializeModels()){
		ConsoleColorizer::PrintRed("Model initializiation failed. Exiting...\n");
		return;
	}

	if(kinX_ == 2){
		TString outline = Form("Beamspot profile is: %s\nBeamspot profile xPar = %f\tyPar = %f\nBeamspot xOffset = %f\tyOffset = %f\n",profile_->ToString().Data(),profile_->GetParX(),profile_->GetParY(),beamspot_->GetXOffset(),beamspot_->GetYOffset());
		ConsoleColorizer::PrintGreen(outline.Data());
		// std::cout << "Beamspot profile is: " << profile_->ToString() << std::endl;
		// std::cout << "Beamspot profile xPar = " << profile_->GetParX() << "\typar = " << profile_->GetParY() << std::endl;
		// std::cout << "Beamspot xOffset = " << beamspot_->GetXOffset() << "\tyOffset = " << beamspot_->GetYOffset() << "\n" << std::endl;
		Simulate2body(infile, outfile);
	} else if(kinX_ == 3){
		TString outline = Form("Beamspot profile is: %s\nBeamspot profile xPar = %f\tyPar = %f\nBeamspot xOffset = %f\tyOffset = %f\n",profile_->ToString().Data(),profile_->GetParX(),profile_->GetParY(),beamspot_->GetXOffset(),beamspot_->GetYOffset());
		ConsoleColorizer::PrintGreen(outline.Data());
		// std::cout << "Beamspot profile is: " << profile_->ToString() << std::endl;
		// std::cout << "Beamspot xPar = " << profile_->GetParX() << "\typar = " << profile_->GetParY() << std::endl;
		// std::cout << "Beamspot xOffset = " << beamspot_->GetXOffset() << "\tyOffset = " << beamspot_->GetYOffset() << "\n" << std::endl;
		Simulate3body(infile,outfile);
	} else if(kinX_ == 4){
		TString outline = Form("Beamspot profile is: %s\nBeamspot profile xPar = %f\tyPar = %f\nBeamspot xOffset = %f\tyOffset = %f\n",profile_->ToString().Data(),profile_->GetParX(),profile_->GetParY(),beamspot_->GetXOffset(),beamspot_->GetYOffset());
		ConsoleColorizer::PrintGreen(outline.Data());
		// std::cout << "Beamspot profile is: " << profile_->ToString() << std::endl;
		// std::cout << "Beamspot xPar = " << profile_->GetParX() << "\typar = " << profile_->GetParY() << std::endl;
		// std::cout << "Beamspot xOffset = " << beamspot_->GetXOffset() << "\tyOffset = " << beamspot_->GetYOffset() << "\n" << std::endl;
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

	//set current targetLoss and deadLayerLoss pointers here! (eventually, this should be done by a config file to avoid remaking every time! In the mean time, we just print out so it is obvious what we are using each time we run!)
	
	//										7Li(3He,4He)6Li
	targetLoss_par1_ = targetLoss_alpha_in_LiF_; //can be set to none, but ejectile is an alpha
	targetLoss_par2_ = targetLoss_6Li_in_LiF_;
	targetLoss_par3_ = targetLoss_none_;
	targetLoss_par4_ = targetLoss_none_;

	deadLayerLoss_par1_ = deadLayerLoss_alpha_; //can be set to none, but ejectile is an alpha
	deadLayerLoss_par2_ = deadLayerLoss_6Li_;
	deadLayerLoss_par3_ = deadLayerLoss_none_;
	deadLayerLoss_par4_ = deadLayerLoss_none_;
	//------------------------------------------------------------------------------------------------

	det2mc det2mcProcessor(SABRE_Array_,
						   SABREARRAY_EnergyResolutionModels_,
						   targetLoss_par1_, targetLoss_par2_,
						   deadLayerLoss_par1_, deadLayerLoss_par2_,
						   beamspot_);

	TString outline = Form("\nPar 1 Target Loss = %s\nPar 2 Target Loss = %s\n\nPar 1 Dead Layer Loss = %s\nPar 2 Dead Layer Loss = %s\n\n",det2mcProcessor.GetToString_TargetLoss1().data(),det2mcProcessor.GetToString_TargetLoss2().data(),det2mcProcessor.GetToString_DeadLayerLoss1().data(),det2mcProcessor.GetToString_DeadLayerLoss2().data());
	ConsoleColorizer::PrintGreen(outline.Data());


	//det2mc det2mcProcessor(SABRE_Array_, SABREARRAY_EnergyResolutionModels_, targetLoss_, deadLayerLoss_, beamspot_);


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

	//set current targetLoss and deadLayerLoss pointers here! (eventually, this should be done by a config file to avoid remaking every time! In the mean time, we just print out so it is obvious what we are using each time we run!)

	//							7Li(3He,4He)6Li (3+)	6Li-> 4He + d
	targetLoss_par1_ = targetLoss_alpha_in_LiF_; 	//can be set to none, but ejectile is an alpha
	targetLoss_par2_ = targetLoss_alpha_in_LiF_; 	//breakup1 is alpha
	targetLoss_par3_ = targetLoss_deuteron_in_LiF_; //breakup2 is deuteron
	targetLoss_par4_ = targetLoss_none_; 			//no particle 4

	deadLayerLoss_par1_ = deadLayerLoss_alpha_; 		//can be set to none, but ejectile is an alpha
	deadLayerLoss_par2_ = deadLayerLoss_alpha_;		//breakup1 is alpha
	deadLayerLoss_par3_ = deadLayerLoss_deuteron_;		//breakup2 is deuteron
	deadLayerLoss_par4_ = deadLayerLoss_none_;			//no particle 4
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


	det3mc det3mcProcessor(SABRE_Array_,
						   SABREARRAY_EnergyResolutionModels_,
						   targetLoss_par1_, targetLoss_par2_, targetLoss_par3_, targetLoss_par4_,
						   deadLayerLoss_par1_, deadLayerLoss_par2_, deadLayerLoss_par3_, deadLayerLoss_par4_,
						   beamspot_);

	//det3mc det3mcProcessor(SABRE_Array_, SABREARRAY_EnergyResolutionModels_, targetLoss_, deadLayerLoss_, beamspot_);
	TString outline = Form("\nPar 1 Target Loss = %s\nPar 2 Target Loss = %s\n\nPar 3 Target Loss = %s\nPar 4 Target Loss = %s\n\nPar 1 Dead Layer Loss = %s\nPar 2 Dead Layer Loss = %s\nPar 3 Dead Layer Loss = %s\nPar 4 Dead Layer Loss = %s\n\n",det3mcProcessor.GetToString_TargetLoss1().data(),det3mcProcessor.GetToString_TargetLoss2().data(),det3mcProcessor.GetToString_TargetLoss3().data(),det3mcProcessor.GetToString_TargetLoss4().data(),det3mcProcessor.GetToString_DeadLayerLoss1().data(),det3mcProcessor.GetToString_DeadLayerLoss2().data(),det3mcProcessor.GetToString_DeadLayerLoss3().data(),det3mcProcessor.GetToString_DeadLayerLoss4().data());
	ConsoleColorizer::PrintGreen(outline.Data());

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
	
	if(!infile.is_open()){
		ConsoleColorizer::PrintRed("Error: Cannot open input file!\n");
		return;
	}

	if(!outfile.is_open()){
		ConsoleColorizer::PrintRed("Error: Cannot open output file!\n");
		return;
	}

	//set current targetLoss and deadLayerLoss pointers here! (eventually, this should be done by a config file to avoid remaking every time! In the mean time, we just print out so it is obvious what we are using each time we run!)

	//							10B(3He,4He)9B		9B -> p + 8Be,	8Be -> 4He + 4He
	targetLoss_par1_ = targetLoss_alpha_in_10B_; 		//can be set to none, but ejectile is an alpha
	targetLoss_par2_ = targetLoss_proton_in_10B_; 		//breakup1 is proton
	targetLoss_par3_ = targetLoss_alpha_in_10B_; 		//breakup2 is alpha
	targetLoss_par4_ = targetLoss_alpha_in_10B_;		//breakup3/final daughter is alpha

	deadLayerLoss_par1_ = deadLayerLoss_alpha_; 			//can be set to none, but ejectile is an alpha
	deadLayerLoss_par2_ = deadLayerLoss_proton_;			//breakup1 is proton
	deadLayerLoss_par3_ = deadLayerLoss_alpha_;			//breakup2 is alpha
	deadLayerLoss_par4_ = deadLayerLoss_alpha_;			//breakup3/final daughter is alpha
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	//							10B(3He,4He)9B		9B -> 4He + 5Li,	5Li -> p + 4He
	// targetLoss_par1_ = targetLoss_alpha_in_10B_; 	//can be set to none, but ejectile is an alpha
	// targetLoss_par2_ = targetLoss_alpha_in_10B_; 	//breakup1 is alpha
	// targetLoss_par3_ = targetLoss_proton_in_10B_; 	//breakup2 is proton
	// targetLoss_par4_ = targetLoss_alpha_in_10B_; 	//breakup3/final daughter is alpha

	// targetLoss_par1_ = deadLayerLoss_alpha_; 		//can be set to none, but ejectile is an alpha
	// targetLoss_par2_ = deadLayerLoss_proton_;		//breakup1 is alpha
	// targetLoss_par3_ = deadLayerLoss_alpha_;			//breakup2 is proton
	// targetLoss_par4_ = deadLayerLoss_alpha_;			//breakup3/final daughter is alpha
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	det4mc det4mcProcessor(SABRE_Array_,
						   SABREARRAY_EnergyResolutionModels_,
						   targetLoss_par1_, targetLoss_par2_, targetLoss_par3_, targetLoss_par4_,
						   deadLayerLoss_par1_, deadLayerLoss_par2_, deadLayerLoss_par3_, deadLayerLoss_par4_,
						   beamspot_);

	TString outline = Form("\nPar 1 Target Loss = %s\nPar 2 Target Loss = %s\n\nPar 3 Target Loss = %s\nPar 4 Target Loss = %s\n\nPar 1 Dead Layer Loss = %s\nPar 2 Dead Layer Loss = %s\nPar 3 Dead Layer Loss = %s\nPar 4 Dead Layer Loss = %s\n\n",det4mcProcessor.GetToString_TargetLoss1().data(),det4mcProcessor.GetToString_TargetLoss2().data(),det4mcProcessor.GetToString_TargetLoss3().data(),det4mcProcessor.GetToString_TargetLoss4().data(),det4mcProcessor.GetToString_DeadLayerLoss1().data(),det4mcProcessor.GetToString_DeadLayerLoss2().data(),det4mcProcessor.GetToString_DeadLayerLoss3().data(),det4mcProcessor.GetToString_DeadLayerLoss4().data());
	ConsoleColorizer::PrintGreen(outline.Data());

	det4mcProcessor.Run(infile,outfile);

	nevents_ = det4mcProcessor.GetNumEvents();
	detectorHits_ = det4mcProcessor.GetDetectorHits();
	hitej_ = det4mcProcessor.GetHitEj();
	hit1_ = det4mcProcessor.GetHit1();
	hit2_ = det4mcProcessor.GetHit2();
	hit3_ = det4mcProcessor.GetHit3();
	hit23_ = det4mcProcessor.GetHitBoth23();
	hitOnlyEj_ = det4mcProcessor.GetHitOnlyEj();
	hitOnly1_ = det4mcProcessor.GetHitOnly1();
	hitOnly2_ = det4mcProcessor.GetHitOnly2();
	hitOnly3_ = det4mcProcessor.GetHitOnly3();
	hitOnly12_ = det4mcProcessor.GetHitOnly12();
	hitOnly23_ = det4mcProcessor.GetHitOnly23();
	hitOnly13_ = det4mcProcessor.GetHitOnly13();
	hitOnly123_ = det4mcProcessor.GetHitOnly123();
	onePartHits_ = det4mcProcessor.GetOnePartHits();
	twoPartHits_ = det4mcProcessor.GetTwoPartHits();
	threePartHits_ = det4mcProcessor.GetThreePartHits();
	fourPartHits_ = det4mcProcessor.GetFourPartHits();
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
		std::cout << "Only bu1: " << hitOnly1_ << "\n";
		std::cout << "Only bu2: " << hitOnly2_ << "\n";
		std::cout << "Only bu3: " << hitOnly3_ << "\n";
		std::cout << "Only bu1 & bu2: " << hitOnly12_ << "\n";
		std::cout << "Only bu2 & bu3: " << hitOnly23_ << "\n";
		std::cout << "Only bu1 & bu3: " << hitOnly13_ << "\n";
		std::cout << "Only bu1, bu2 & bu3: " << hitOnly123_ << std::endl;
	}
	for (size_t i = 0; i < detectorHits_.size(); i++) {
		std::cout << "Detector_" << i << " total hits = " << detectorHits_[i] << "\n";
	}
	std::cout << "Output file: " << detOutputFilename_ << "\n\n";
}
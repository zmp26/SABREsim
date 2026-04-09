#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <vector>

std::string MakeOutputName(const char* infile, const std::string& suffix){
	std::string name(infile);
	size_t pos = name.rfind(".root");
	if(pos != std::string::npos){
		name.insert(pos,suffix);
	} else {
		name += suffix + ".root";
	}
	return name;
}

void SABREslicer(const char* infile, const char* intree = "SABREsim", int omniscienceSort=1){
	/*
	
		SABREslicer arguments:

		const char* infile 					path to ROOT file containing "SABREsim" TTree output from SABREsim/bin/SABREsim
		const char* intree					name of TTree in ROOT file. Default value is "SABREsim" for SABREsim output
		int omniscienceSort					int to decide whether to sort particle arrays by particle ID (1), randomize order(2), or keep same(0)



		A note on omniscienceSort:

			SABREsim det2/3/4mc classes currently iterate through the SABRE_Array object from SABRE0->SABRE4 and checks *all* particles
			against that detector before moving on to the next detector in the array. As a consequence of this behavior, the output of
			SABREsim is stored in a potentially different order when compared to the input of the kin2/3/4mc file. This means that the
			results will be biased by the laboratory angles of the simulated break-up particles, with those more likely to be in SABRE0
			to be near the front of the list and the others near the back of the list. If nothing is done about this, then our cases in
			SABREsim/analyze/SABREanalyze_multX.cxx codes will have the distribution of events across the pre-defined cases biased in
			the same fashion. There are two ways I can think of to get around this:

				1) Sort the particles by particleID before saving to the output multiplicity TTrees. This preserves omniscience and lets
				   us understand our cases a lot better for simulated data. Use omniscienceSort=1 for this option. (DEFAULT)

				2) Randomize the particle orders before saving to the output multiplicity TTrees. This will prevent the aforementioned bias
				   and *could* be more representative of the data since the experimental trigger was the SPS which just opened a coincidence
				   window to record SABRE detections regardless of which specific detector they hit. Use omniscienceSort=2 for this option.

			We may also keep the order as output by SABREsim. Use omniscienceSort=0 for this option, although any non-1 and non-2 integer works.
			

	*/
	TFile *fin = new TFile(infile, "READ");
	if(!fin || fin->IsZombie()){
		std::cerr << "Error opening " << infile << std::endl;
		return;
	}

	TTree *tin = (TTree*)fin->Get(intree);
	if(!tin){
		std::cerr << "TTree " << intree << " not found in " << infile << std::endl;
		return;
	}

	std::string out1 = MakeOutputName(infile, "_mult1");
	std::string out2 = MakeOutputName(infile, "_mult2");
	std::string out3 = MakeOutputName(infile, "_mult3");

	TFile *f1 = new TFile(out1.c_str(), "RECREATE");
	TTree *t1 = new TTree("mult1","mult1");

	TFile *f2 = new TFile(out2.c_str(), "RECREATE");
	TTree *t2 = new TTree("mult2", "mult2");

	TFile *f3 = new TFile(out3.c_str(), "RECREATE");
	TTree *t3 = new TTree("mult3", "mult3");

	//branch variables for new trees:
	int eventnum;
	double ExEplaceholder, SPSEnergy, SPSTheta, SPSPhi;
	double SR1, SW1, SRE1, SWE1, Stheta1, Sphi1;
	double SR2, SW2, SRE2, SWE2, Stheta2, Sphi2;
	double SR3, SW3, SRE3, SWE3, Stheta3, Sphi3;

	double* const SR[3]        = { &SR1, &SR2, &SR3 };
	double* const SW[3]        = { &SW1, &SW2, &SW3 };
	double* const SRE[3]       = { &SRE1, &SRE2, &SRE3 };
	double* const SWE[3]       = { &SWE1, &SWE2, &SWE3 };
	double* const Stheta[3]    = { &Stheta1, &Stheta2, &Stheta3 };
	double* const Sphi[3]      = { &Sphi1, &Sphi2, &Sphi3 };

	t1->Branch("eventnum", &eventnum, "eventnum/I");
	t1->Branch("ExEPlaceHolder", &ExEplaceholder, "ExEPlaceHolder/D");
	t1->Branch("SPSEnergy", &SPSEnergy, "SPSEnergy/D");
	t1->Branch("SPSTheta", &SPSTheta, "SPSTheta/D");
	t1->Branch("SPSPhi", &SPSPhi, "SPSPhi/D");
	t1->Branch("SabreRing_hit1", &SR1, "SabreRing_hit1/D");
	t1->Branch("SabreWedge_hit1", &SW1, "SabreWedge_hit1/D");
	t1->Branch("SabreRingEnergy_hit1", &SRE1, "SabreRingEnergy_hit1/D");
	t1->Branch("SabreWedgeEnergy_hit1", &SWE1, "SabreWedgeEnergy_hit1/D");
	t1->Branch("thetalab_hit1", &Stheta1, "thetalab_hit1/D");
	t1->Branch("philab_hit1", &Sphi1, "philab_hit1/D");


	t2->Branch("eventnum", &eventnum, "eventnum/I");
	t2->Branch("ExEPlaceHolder", &ExEplaceholder, "ExEPlaceHolder/D");
	t2->Branch("SPSEnergy", &SPSEnergy, "SPSEnergy/D");
	t2->Branch("SPSTheta", &SPSTheta, "SPSTheta/D");
	t2->Branch("SPSPhi", &SPSPhi, "SPSPhi/D");
	t2->Branch("SabreRing_hit1", &SR1, "SabreRing_hit1/D");
	t2->Branch("SabreWedge_hit1", &SW1, "SabreWedge_hit1/D");
	t2->Branch("SabreRingEnergy_hit1", &SRE1, "SabreRingEnergy_hit1/D");
	t2->Branch("SabreWedgeEnergy_hit1", &SWE1, "SabreWedgeEnergy_hit1/D");
	t2->Branch("thetalab_hit1", &Stheta1, "thetalab_hit1/D");
	t2->Branch("philab_hit1", &Sphi1, "philab_hit1/D");
	t2->Branch("SabreRing_hit2", &SR2, "SabreRing_hit2/D");
	t2->Branch("SabreWedge_hit2", &SW2, "SabreWedge_hit2/D");
	t2->Branch("SabreRingEnergy_hit2", &SRE2, "SabreRingEnergy_hit2/D");
	t2->Branch("SabreWedgeEnergy_hit2", &SWE2, "SabreWedgeEnergy_hit2/D");
	t2->Branch("thetalab_hit2", &Stheta2, "thetalab_hit2/D");
	t2->Branch("philab_hit2", &Sphi2, "philab_hit2/D");


	t3->Branch("eventnum", &eventnum, "eventnum/I");
	t3->Branch("ExEPlaceHolder", &ExEplaceholder, "ExEPlaceHolder/D");
	t3->Branch("SPSEnergy", &SPSEnergy, "SPSEnergy/D");
	t3->Branch("SPSTheta", &SPSTheta, "SPSTheta/D");
	t3->Branch("SPSPhi", &SPSPhi, "SPSPhi/D");
	t3->Branch("SabreRing_hit1", &SR1, "SabreRing_hit1/D");
	t3->Branch("SabreWedge_hit1", &SW1, "SabreWedge_hit1/D");
	t3->Branch("SabreRingEnergy_hit1", &SRE1, "SabreRingEnergy_hit1/D");
	t3->Branch("SabreWedgeEnergy_hit1", &SWE1, "SabreWedgeEnergy_hit1/D");
	t3->Branch("thetalab_hit1", &Stheta1, "thetalab_hit1/D");
	t3->Branch("philab_hit1", &Sphi1, "philab_hit1/D");
	t3->Branch("SabreRing_hit2", &SR2, "SabreRing_hit2/D");
	t3->Branch("SabreWedge_hit2", &SW2, "SabreWedge_hit2/D");
	t3->Branch("SabreRingEnergy_hit2", &SRE2, "SabreRingEnergy_hit2/D");
	t3->Branch("SabreWedgeEnergy_hit2", &SWE2, "SabreWedgeEnergy_hit2/D");
	t3->Branch("thetalab_hit2", &Stheta2, "thetalab_hit2/D");
	t3->Branch("philab_hit2", &Sphi2, "philab_hit2/D");
	t3->Branch("SabreRing_hit3", &SR3, "SabreRing_hit3/D");
	t3->Branch("SabreWedge_hit3", &SW3, "SabreWedge_hit3/D");
	t3->Branch("SabreRingEnergy_hit3", &SRE3, "SabreRingEnergy_hit3/D");
	t3->Branch("SabreWedgeEnergy_hit3", &SWE3, "SabreWedgeEnergy_hit3/D");
	t3->Branch("thetalab_hit3", &Stheta3, "thetalab_hit3/D");
	t3->Branch("philab_hit3", &Sphi3, "philab_hit3/D");

	//prepare variables to read value of tin entry
	double tin_kine[4], tin_kintheta[4], tin_kinphi[4];
	double tin_reactionOrigin[3];
	int tin_multiplicity;
	int tin_particleID[4], tin_detectorID[4], tin_ringChannel[4], tin_wedgeChannel[4], tin_localRing[4], tin_localWedge[4];
	double tin_ringEnergy[4], tin_wedgeEnergy[4], tin_ringTheta[4], tin_wedgePhi[4], tin_localx[4], tin_localy[4];

	tin->SetBranchAddress("kin_e", tin_kine);
	tin->SetBranchAddress("kin_theta", tin_kintheta);
	tin->SetBranchAddress("kin_phi", tin_kinphi);
	tin->SetBranchAddress("reactionOrigin", tin_reactionOrigin);
	tin->SetBranchAddress("numHits", &tin_multiplicity);
	tin->SetBranchAddress("particleID", tin_particleID);
	tin->SetBranchAddress("detectorID", tin_detectorID);
	tin->SetBranchAddress("ringChannel", tin_ringChannel);
	tin->SetBranchAddress("wedgeChannel", tin_wedgeChannel);
	tin->SetBranchAddress("localRing", tin_localRing);
	tin->SetBranchAddress("localWedge", tin_localWedge);
	tin->SetBranchAddress("ringEnergy", tin_ringEnergy);
	tin->SetBranchAddress("wedgeEnergy", tin_wedgeEnergy);
	tin->SetBranchAddress("ringTheta", tin_ringTheta);
	tin->SetBranchAddress("wedgePhi", tin_wedgePhi);
	tin->SetBranchAddress("localx", tin_localx);
	tin->SetBranchAddress("localy", tin_localy);

	//random number generator for shuffling
	std::random_device rd;
	std::mt19937 g(rd());

	//loop over input tree:
	Long64_t nentries = tin->GetEntries();
	for(Long64_t i=0; i<nentries; i++){
	
		tin->GetEntry(i);

		eventnum = (int) i;
		ExEplaceholder = -666.666; //this is to match output format of Rachel2Root for experimental data, but is not used. Vestigial branch.

		for(int k=0; k<3; k++){
			*SR[k] = -666.;
			*SW[k] = -666.;
			*SRE[k] = -666.;
			*SWE[k] = -666.;
			*Stheta[k] = -666.;
			*Sphi[k] = -666.;
		}

		SPSEnergy = tin_kine[0];
		SPSTheta = tin_kintheta[0];
		SPSPhi = tin_kinphi[0];

		//index list for shuffling/sorting:
		int nParticles = std::min(tin_multiplicity,3);
		std::vector<int> indices(nParticles);
		for(int j=0; j<nParticles; j++) indices[j] = j;

		//apply sorting randomization here:
		if(omniscienceSort == 1){
			std::sort(indices.begin(), indices.end(), [&](int a, int b){
				return tin_particleID[a] < tin_particleID[b];
			});
		} else if(omniscienceSort == 2){
			std::shuffle(indices.begin(), indices.end(), g);
		}

		for(int k=0; k<nParticles; k++){
			int j = indices[k];
			*SR[k] = tin_localRing[j];
			*SW[k] = tin_localWedge[j];
			*SRE[k] = tin_ringEnergy[j];
			*SWE[k] = tin_wedgeEnergy[j];
			*Stheta[k] = tin_ringTheta[j];
			*Sphi[k] = tin_wedgePhi[j];
		}

		if(nParticles == 1){//SABRE multiplicity 1

			//simply fill the tree if the multiplicity is correct:
			t1->Fill();

		} else if(nParticles == 2){//SABRE multiplicity 2

			//simply fill the tree if the multiplicity is correct:
			t2->Fill();

		} else if(nParticles == 3){//SABRE multiplicity 3

			//simply fill the tree if the multiplicity is correct:
			t3->Fill();
		}

	}

	//write the trees to their respective files here:
	f1->cd();
	t1->Write();
	f1->Close();

	f2->cd();
	t2->Write();
	f2->Close();

	f3->cd();
	t3->Write();
	f3->Close();

	std::cout << "\nCreated:\n" 
			  << "\t" << out1 << "\n"
			  << "\t" << out2 << "\n"
			  << "\t" << out3 << "\n\n";

	fin->Close();

	delete f1;
	delete f2;
	delete f3;
	delete fin;
}



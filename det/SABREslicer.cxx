#include "TFile.h"
#include "TTree.h"
#include <iostream>

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

void SABREslicer(const char* infile, const char* intree = "SABREsim"){
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

	int numHits;
	tin->SetBranchAddress("numHits", &numHits);

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

	//loop over input tree:
	Long64_t nentries = tin->GetEntries();
	for(Long64_t i=0; i<nentries; i++){
	
		tin->GetEntry(i);

		eventnum = (int) i;
		ExEplaceholder = -666.666;

		SPSEnergy = tin_kine[0];
		SPSTheta = tin_kintheta[0];
		SPSPhi = tin_kinphi[0];

		SR1 = tin_localRing[0];
		SW1 = tin_localWedge[0];
		SRE1 = tin_ringEnergy[0];
		SWE1 = tin_wedgeEnergy[0];
		Stheta1 = tin_ringTheta[0];
		Sphi1 = tin_wedgePhi[0];

		SR2 = tin_localRing[1];
		SW2 = tin_localWedge[1];
		SRE2 = tin_ringEnergy[1];
		SWE2 = tin_wedgeEnergy[1];
		Stheta2 = tin_ringTheta[1];
		Sphi2 = tin_wedgePhi[1];

		SR3 = tin_localRing[2];
		SW3 = tin_localWedge[2];
		SRE3 = tin_ringEnergy[2];
		SWE3 = tin_wedgeEnergy[2];
		Stheta3 = tin_ringTheta[2];
		Sphi3 = tin_wedgePhi[2];

		if(tin_multiplicity == 1){//SABRE multiplicity 1

			//simply fill the tree if the multiplicity is correct:
			t1->Fill();

		} else if(tin_multiplicity == 2){//SABRE multiplicity 2

			//simply fill the tree if the multiplicity is correct:
			t2->Fill();

		} else if(tin_multiplicity == 3){//SABRE multiplicity 3

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
}



#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <cstring>

std::string MakePairOutputName(const char* infile, const std::string& pairSuffix) {
    std::string name(infile);
    size_t pos = name.rfind(".root");
    if (pos != std::string::npos) {
        name.insert(pos, "_" + pairSuffix);
    } else {
        name += "_" + pairSuffix + ".root";
    }
    return name;
}

void SABRE_ParticleID_Splitter(const char* infile, const char* intree = "SABREsim") {
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error opening " << infile << std::endl;
        return;
    }

    TTree *tin = (TTree*)fin->Get(intree);
    if (!tin) {
        std::cerr << "TTree " << intree << " not found!" << std::endl;
        fin->Close();
        return;
    }

    // 3 Unique files representing the pure, isolated combinations
    std::string outDE_name = MakePairOutputName(infile, "mult2_de"); // Missing c (particle 2)
    std::string outCD_name = MakePairOutputName(infile, "mult2_cd"); // Missing e (particle 4)
    std::string outCE_name = MakePairOutputName(infile, "mult2_ce"); // Missing d (particle 3)

    TFile *f_de = new TFile(outDE_name.c_str(), "RECREATE");
    TTree *t_de = new TTree("mult2", "Detected d+e (Missing c)");

    TFile *f_cd = new TFile(outCD_name.c_str(), "RECREATE");
    TTree *t_cd = new TTree("mult2", "Detected c+d (Missing e)");

    TFile *f_ce = new TFile(outCE_name.c_str(), "RECREATE");
    TTree *t_ce = new TTree("mult2", "Detected c+e (Missing d)");

    // Branch Variables for flat tree output
    int eventnum;
    double kine[4], kintheta[4], kinphi[4], reactionOrigin[3];
    double ExE;
    bool EjInSPS;
    double SPSEnergy, SPSTheta, SPSPhi;

    int SR1, SW1; double SRE1, SWE1, Stheta1, Sphi1;
    int SR2, SW2; double SRE2, SWE2, Stheta2, Sphi2;
    int key[2];

    auto BindBranches = [&](TTree* t) {
        t->Branch("eventnum", &eventnum, "eventnum/I");
        t->Branch("kin_e", kine, "kin_e[4]/D");
        t->Branch("kin_theta", kintheta, "kin_theta[4]/D");
        t->Branch("kin_phi", kinphi, "kin_phi[4]/D");
        t->Branch("reactionOrigin", reactionOrigin, "reactionOrigin[3]/D");
        t->Branch("ExE", &ExE, "ExE/D");
        t->Branch("EjInSPS", &EjInSPS, "EjInSPS/O");
        t->Branch("SPSEnergy", &SPSEnergy, "SPSEnergy/D");
        t->Branch("SPSTheta", &SPSTheta, "SPSTheta/D");
        t->Branch("SPSPhi", &SPSPhi, "SPSPhi/D");

        t->Branch("SabreRing_hit1", &SR1, "SabreRing_hit1/I");
        t->Branch("SabreWedge_hit1", &SW1, "SabreWedge_hit1/I");
        t->Branch("SabreRingEnergy_hit1", &SRE1, "SabreRingEnergy_hit1/D");
        t->Branch("SabreWedgeEnergy_hit1", &SWE1, "SabreWedgeEnergy_hit1/D");
        t->Branch("thetalab_hit1", &Stheta1, "thetalab_hit1/D");
        t->Branch("philab_hit1", &Sphi1, "philab_hit1/D");
        t->Branch("key_hit1", &key[0], "key_hit1/I");

        t->Branch("SabreRing_hit2", &SR2, "SabreRing_hit2/I");
        t->Branch("SabreWedge_hit2", &SW2, "SabreWedge_hit2/I");
        t->Branch("SabreRingEnergy_hit2", &SRE2, "SabreRingEnergy_hit2/D");
        t->Branch("SabreWedgeEnergy_hit2", &SWE2, "SabreWedgeEnergy_hit2/D");
        t->Branch("thetalab_hit2", &Stheta2, "thetalab_hit2/D");
        t->Branch("philab_hit2", &Sphi2, "philab_hit2/D");
        t->Branch("key_hit2", &key[1], "key_hit2/I");
    };

    BindBranches(t_de); BindBranches(t_cd); BindBranches(t_ce);

    // Tree extraction variables
    double tin_kine[4], tin_kintheta[4], tin_kinphi[4], tin_reactionOrigin[3];
    bool tin_EjInSPS;
    double tin_SPSEnergy, tin_SPSTheta, tin_SPSPhi, tin_ExE;
    int tin_multiplicity;
    int tin_particleID[4], tin_localRing[4], tin_localWedge[4];
    double tin_ringEnergy[4], tin_wedgeEnergy[4], tin_ringTheta[4], tin_wedgePhi[4];

    tin->SetBranchAddress("kin_e", tin_kine);
    tin->SetBranchAddress("kin_theta", tin_kintheta);
    tin->SetBranchAddress("kin_phi", tin_kinphi);
    tin->SetBranchAddress("reactionOrigin", tin_reactionOrigin);
    tin->SetBranchAddress("EjInSPS", &tin_EjInSPS);
    tin->SetBranchAddress("SPSEnergy", &tin_SPSEnergy);
    tin->SetBranchAddress("SPSTheta", &tin_SPSTheta);
    tin->SetBranchAddress("SPSPhi", &tin_SPSPhi);
    tin->SetBranchAddress("ExE", &tin_ExE);
    tin->SetBranchAddress("numHits", &tin_multiplicity);
    tin->SetBranchAddress("particleID", tin_particleID);
    tin->SetBranchAddress("localRing", tin_localRing);
    tin->SetBranchAddress("localWedge", tin_localWedge);
    tin->SetBranchAddress("ringEnergy", tin_ringEnergy);
    tin->SetBranchAddress("wedgeEnergy", tin_wedgeEnergy);
    tin->SetBranchAddress("ringTheta", tin_ringTheta);
    tin->SetBranchAddress("wedgePhi", tin_wedgePhi);

    Long64_t nentries = tin->GetEntries();
    
    for (Long64_t i = 0; i < nentries; i++) {
        tin->GetEntry(i);

        if (tin_multiplicity != 2) continue;

        eventnum = (int)i;
        ExE = tin_ExE;
        EjInSPS = tin_EjInSPS;
        SPSEnergy = tin_SPSEnergy;
        SPSTheta = tin_SPSTheta;
        SPSPhi = tin_SPSPhi;

        std::memcpy(kine, tin_kine, sizeof(kine));
        std::memcpy(kintheta, tin_kintheta, sizeof(kintheta));
        std::memcpy(kinphi, tin_kinphi, sizeof(kinphi));
        std::memcpy(reactionOrigin, tin_reactionOrigin, sizeof(reactionOrigin));

        int id1 = tin_particleID[0]; 
        int id2 = tin_particleID[1]; 

        // Retain native ordering index for downstream mapping consistency
        key[0] = 0;
        key[1] = 1;

        SR1 = tin_localRing[0];     SW1 = tin_localWedge[0];
        SRE1 = tin_ringEnergy[0];   SWE1 = tin_wedgeEnergy[0];
        Stheta1 = tin_ringTheta[0]; Sphi1 = tin_wedgePhi[0];

        SR2 = tin_localRing[1];     SW2 = tin_localWedge[1];
        SRE2 = tin_ringEnergy[1];   SWE2 = tin_wedgeEnergy[1];
        Stheta2 = tin_ringTheta[1]; Sphi2 = tin_wedgePhi[1];

        // Direct sorting based on the natural loop output sequence of det4mc
        if (id1 == 3 && id2 == 4) {
            t_de->Fill();
        } else if (id1 == 2 && id2 == 3) {
            t_cd->Fill();
        } else if (id1 == 2 && id2 == 4) {
            t_ce->Fill();
        }
    }

    f_de->cd(); t_de->Write(); f_de->Close(); delete f_de;
    f_cd->cd(); t_cd->Write(); f_cd->Close(); delete f_cd;
    f_ce->cd(); t_ce->Write(); f_ce->Close(); delete f_ce;
    
    fin->Close();
    delete fin;

    std::cout << "Splitting complete. Subfiles maintain native simulation ordering (ascending particle ID numbers)." << std::endl;
}
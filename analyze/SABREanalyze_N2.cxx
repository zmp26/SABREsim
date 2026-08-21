
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <chrono>
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "permHisto_mult2.h"
#include "permHisto_mult2.cpp"
#include "SABREPID_N2_M2.h"
#include "SABREPID_N2_M2.cpp"

void Li7ha_SABREPID_N2_M2(const char* input_filename){
	std::string s = input_filename;
	size_t last_dot = s.find_last_of(".");
	std::string stem = (last_dot == std::string::npos ? s : s.substr(0, last_dot));
	std::string out_str = stem + "_SABREPID_N2_M2.root";

	TMassTable fMassTable;
	fMassTable.Init("/home/zachpurcell/masstable/masstable.dat");

	double massgs_6Li = fMassTable.GetNuclearMassMeV("Li",6);
	double massgs_4He = fMassTable.GetNuclearMassMeV("He",4);
	double massgs_2H = fMassTable.GetNuclearMassMeV("H",2);

	PIDHypothesis_N2 hypothesis;
	hypothesis.mass_target = fMassTable.GetNuclearMassMeV("Li",7);
	hypothesis.mass_beam = fMassTable.GetNuclearMassMeV("He",3);
	hypothesis.mass_ejectile = massgs_4He;
	hypothesis.beamEnergyMeV = 7.5;

	hypothesis.name = "ad";
	hypothesis.final_masses[0] = massgs_4He;
	hypothesis.final_particles[0] = "#alpha";
	hypothesis.final_masses[1] = massgs_2H;
	hypothesis.final_particles[1] = "d";

	TFile *infile = TFile::Open(input_filename, "READ");
	if(!infile || infile->IsZombie()){
		std::cerr << "Error: Cannot open input file " << input_filename << std::endl;
		return;
	}

	TTree *intree = (TTree*)infile->Get("mult2");
	if(!intree){
		std::cerr << "Error: cannot get TTree 'mult2' from " << input_filename << std::endl;
		infile->Close();
		return;
	}

	double Ex, SPSE, SPSTheta, SPSPhi;
	intree->SetBranchAddress("ExE", &Ex);
	intree->SetBranchAddress("SPSEnergy", &SPSE);
	intree->SetBranchAddress("SPSTheta", &SPSTheta);
	intree->SetBranchAddress("SPSPhi", &SPSPhi);

	double E[3], theta[3], phi[3];
	intree->SetBranchAddress("SabreRingEnergy_hit1", &E[0]);
	intree->SetBranchAddress("thetalab_hit1", &theta[0]);
	intree->SetBranchAddress("philab_hit1", &phi[0]);

	intree->SetBranchAddress("SabreRingEnergy_hit2", &E[1]);
	intree->SetBranchAddress("thetalab_hit2", &theta[1]);
	intree->SetBranchAddress("philab_hit2", &phi[1]);

	long numentries = intree->GetEntries();

	TFile *outfile = new TFile(out_str.c_str(), "RECREATE");
	//TDirectory *sliceDir = outfile->mkdir("slices");
	outfile->cd();

	//eventually add histograms below:
	TH2D *hRecoilExSPS_vs_RecoilExSABRE = new TH2D("hRecoilExSPS_vs_RecoilExSABRE", "Recoil Ex SPS vs Recoil Ex SABRE;SABRE;SPS", 100, 0, 10, 100, 0, 10);
	hRecoilExSPS_vs_RecoilExSABRE->SetDirectory(outfile);

	TH2D *hRecoilFlipExSPS_vs_RecoilExSABRE = new TH2D("hRecoilFlipExSPS_vs_RecoilExSABRE", "Recoil Flip Ex SPS vs Recoil Ex SABRE;SABRE;SPS", 100, 0, 10, 100, 0, 10);
	hRecoilFlipExSPS_vs_RecoilExSABRE->SetDirectory(outfile);

	TH1D *hRecoilExSABRE = new TH1D("hRecoilExSABRE", "Recoil EX SABRE", 100, 0, 10);
	hRecoilExSABRE->SetDirectory(outfile);

	TH1D *hRecoilFlipExSABRE = new TH1D("hRecoilFlipExSABRE", "Recoil EX SABRE", 100, 0, 10);
	hRecoilFlipExSABRE->SetDirectory(outfile);

	TH1D *hRecoilBothExSABRE = new TH1D("hRecoilBothExSABRE", "Recoil Ex (both assignments)", 100, 0, 10);
	hRecoilBothExSABRE->SetDirectory(outfile);

	TH1D *hAlphaPicks = new TH1D("hAlphaPicks", "Hit Index Assigned to Alpha", 2, -0.5, 1.5);
	hAlphaPicks->SetDirectory(outfile);

	TH1D *hThetaHAlpha = new TH1D("hThetaHAlpha", "#theta^{h}_{#alpha}", 360, 0., 180.);
	hThetaHAlpha->SetDirectory(outfile);

	TH1D *hThetaHDeuteron = new TH1D("hThetaHDeuteron", "#theta^{h}_{d}", 360, 0., 180.);
	hThetaHDeuteron->SetDirectory(outfile);

	TH1D *hThetaHSum = new TH1D("hThetaHSum", "#theta^{h}_{sum}", 360, 0., 360.);
	hThetaHSum->SetDirectory(outfile);

	TH1D *hCosThetaHAlpha = new TH1D("hCosThetaHAlpha", "cos(#theta^{h}_{#alpha})", 100, -1., 1.);
	hCosThetaHAlpha->SetDirectory(outfile);

	TH1D *hCosThetaHDeuteron = new TH1D("hCosThetaHDeuteron", "cos(#theta^{h}_{d})", 100, -1., 1.);
	hCosThetaHDeuteron->SetDirectory(outfile);

	SABREPID_N2_M2 pidSolver;
	pidSolver.SetHypothesis(hypothesis);
	pidSolver.SetResolution(0.05, 1.0, 1.0);
	pidSolver.SetSPSResolution(0.015, 0.5, 0.5);
	pidSolver.SetChi2Cut(10.);
	pidSolver.InitDiagnostics(outfile);

	TTree *outtree = new TTree("PID_N2_M2", "PID_N2_M2");
	outtree->SetDirectory(outfile);

	int bestPermIndex;
	double bestChi2;
	bool passesCut;
	int alpha_hit_index, deuteron_hit_index;

	TLorentzVector alpha, alphaflip, deuteron, deuteronflip, recoil, recoilflip;
	double residual_px, residual_py, residual_pz, residual_Pmag;

	double ExSPS, recoilEx, recoilFlipEx;
	double cosThetaH_alpha, cosThetaH_deuteron, thetaH_alpha, thetaH_deuteron;

	outtree->Branch("bestPermIndex", &bestPermIndex, "bestPermIndex/I");
	outtree->Branch("bestChi2", &bestChi2, "bestChi2/D");
	outtree->Branch("passesCut", &passesCut, "passesCut/O");

	outtree->Branch("alphaHitIndex", &alpha_hit_index, "alphaHitIndex/I");
	outtree->Branch("deuteronHitIndex", &deuteron_hit_index, "deuteronHitIndex/I");

	outtree->Branch("P4_alpha", &alpha);
	outtree->Branch("P4_deuteron", &deuteron);
	outtree->Branch("P4_recoil", &recoil);

	outtree->Branch("ExSPS", &ExSPS, "ExSPS/D");
	outtree->Branch("RecEx", &recoilEx, "recoilEx/D");

	outtree->Branch("residual_px", &residual_px, "residual_px/D");
	outtree->Branch("residual_py", &residual_py, "residual_py/D");
	outtree->Branch("residual_pz", &residual_pz, "residual_pz/D");
	outtree->Branch("residual_Pmag", &residual_Pmag, "residual_Pmag/D");

	outtree->Branch("thetaH_alpha", &thetaH_alpha, "thetaH_alpha/D");
	outtree->Branch("thetaH_deuteron", &thetaH_deuteron, "thetaH_deuteron/D");
	outtree->Branch("cosThetaH_alpha", &cosThetaH_alpha, "cosThetaH_alpha/D");
	outtree->Branch("cosThetaH_deuteron", &cosThetaH_deuteron, "cosThetaH_deuteron/D");


	for(long i=0; i<numentries; i++){

		intree->GetEntry(i);

		PIDResult_N2_M2 res = pidSolver.EvaluateEvent(E, theta, phi, SPSE, SPSTheta, SPSPhi);

		ExSPS = Ex;

		bestPermIndex = res.bestChi2Index;
		bestChi2 = res.bestChi2;
		passesCut = res.passesCut;

		alpha_hit_index = res.hit_indices[0];
		deuteron_hit_index = res.hit_indices[1];

		residual_px = res.missing_px;
		residual_py = res.missing_py;
		residual_pz = res.missing_pz;
		residual_Pmag = res.missing_Pmag;

		if(bestPermIndex >= 0){
			auto buildP4 = [](double E, double theta, double phi, double m){
				double p = std::sqrt(E*(E+2.*m));
				double rad_th = theta*M_PI/180.;
				double rad_ph = phi*M_PI/180.;
				return TLorentzVector(
						p*std::sin(rad_th)*std::cos(rad_ph),
						p*std::sin(rad_th)*std::sin(rad_ph),
						p*std::cos(rad_th),
						E+m
					);
			};


			//fill pick histograms here:
			hAlphaPicks->Fill(alpha_hit_index);

			alpha = buildP4(E[alpha_hit_index], theta[alpha_hit_index], phi[alpha_hit_index], massgs_4He);
			alphaflip = buildP4(E[deuteron_hit_index], theta[deuteron_hit_index], phi[deuteron_hit_index], massgs_4He);
			deuteron = buildP4(E[deuteron_hit_index], theta[deuteron_hit_index], phi[deuteron_hit_index], massgs_2H);
			deuteronflip = buildP4(E[alpha_hit_index], theta[alpha_hit_index], phi[alpha_hit_index], massgs_2H);
			recoil = alpha + deuteron;
			recoilflip = alphaflip + deuteronflip;
			TVector3 betaLabToCM = -recoil.BoostVector();

			TLorentzVector alphaCM = alpha;
			TLorentzVector deuteronCM = deuteron;
			alphaCM.Boost(betaLabToCM);
			deuteronCM.Boost(betaLabToCM);

			thetaH_alpha = -666.;
			cosThetaH_alpha = -666.;
			thetaH_deuteron = -666.;
			cosThetaH_deuteron = -666.;

			if(betaLabToCM.Mag() > 1e-12 && alphaCM.P() > 1e-12 && deuteronCM.P() > 1e-12){
				TVector3 helicityAxis = betaLabToCM.Unit();

				cosThetaH_alpha = std::max(-1., std::min(1., alphaCM.Vect().Unit().Dot(helicityAxis)));
				cosThetaH_deuteron = std::max(-1., std::min(1., deuteronCM.Vect().Unit().Dot(helicityAxis)));
			
				thetaH_alpha = std::acos(cosThetaH_alpha)*RADDEG;
				thetaH_deuteron = std::acos(cosThetaH_deuteron)*RADDEG;
			}

			hThetaHAlpha->Fill(thetaH_alpha);
			hCosThetaHAlpha->Fill(cosThetaH_alpha);
			hThetaHDeuteron->Fill(thetaH_deuteron);
			hCosThetaHDeuteron->Fill(cosThetaH_deuteron);
			hThetaHSum->Fill(thetaH_alpha + thetaH_deuteron);

			recoilEx = recoil.M() - massgs_6Li;
			recoilFlipEx = recoilflip.M() - massgs_6Li;
			hRecoilExSPS_vs_RecoilExSABRE->Fill(recoilEx, Ex);
			hRecoilFlipExSPS_vs_RecoilExSABRE->Fill(recoilFlipEx, Ex);
			hRecoilExSABRE->Fill(recoilEx);
			hRecoilFlipExSABRE->Fill(recoilFlipEx);
			hRecoilBothExSABRE->Fill(recoilEx); hRecoilBothExSABRE->Fill(recoilFlipEx);

			outtree->Fill();

		}

		if(i % 10000 == 0){
			std::cout << "Processed " << i << " events..." << std::endl;
		}
	}

	infile->Close();
	outfile->cd();
	//hRecoilExSPS_vs_RecoilExSABRE->Write();
	outfile->Write();
	outfile->Close();
	delete outfile;

	std::cout << "Finished! Processed " << numentries << " events. Output stored in " << out_str.c_str() << std::endl;
}
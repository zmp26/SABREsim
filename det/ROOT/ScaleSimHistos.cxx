#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include <iostream>

void scaleHist(const char* inputFileName, const char* inHistName, const char* outHistName, const char* outputFileName, double scaleFactor){

	TFile* infile = TFile::Open(inputFileName, "READ");
	if(!infile || infile->IsZombie()){
		std::cerr << "bad root file\n";
		return;
	}

	TH1D *h = dynamic_cast<TH1D*>(infile->Get(inHistName));
	if(!h){
		std::cerr << "bad histo or histo not found\n";
		infile->Close();
		return;
	}

	TH1D* hScaled = dynamic_cast<TH1D*>(h->Clone(Form("%s_scaled",outHistName)));

	for(int i = 1; i <= hScaled->GetNbinsX(); i++){
		hScaled->SetBinContent(i, hScaled->GetBinContent(i)*scaleFactor);
		//hScaled->SetBinError(i, hScaled->GetBinError(i)*scaleFactor);
	}

	hScaled->GetXaxis()->CenterTitle();
	hScaled->GetYaxis()->CenterTitle();
	hScaled->GetXaxis()->SetTitle("^{6}Li Excitation Energy (MeV)");
	hScaled->GetYaxis()->SetTitle("Counts / 25 keV");
	hScaled->SetTitle("");
	hScaled->SetStats(0);
	hScaled->GetXaxis()->SetRangeUser(0,7);
	hScaled->Draw();

	// TFile* outfile = TFile::Open(outputFileName, "RECREATE");
	// if(!outfile || outfile->IsZombie()){
	// 	std::cerr << "bad out file\n";
	// 	infile->Close();
	// 	return;
	// }

	// outfile->cd();
	// hScaled->Write();

	// outfile->Close();
	// infile->Close();

}

/*


commands for:

	6Li 2.186 MeV a+d:
		
		IMM:	
		scaleHist("/mnt/e/SABREsim/det/kin3mc/kin3mc_7Li3He4He6Li2186keV_4He2H_7500keV_theta178218_phi_-2.125_2.125_detPlothistos.root", "2par/h2par_RecInvMassExE", "hIMM_RecoilExE_SIMSCALED", "/mnt/e/SABREsim/det/kin3mc/kin3mc_7Li3He4He6Li2186keV_4He2H_7500keV_theta178218_phi_-2.125_2.125_histos_IMM_Scaled.root", 6228./1247004.)

		MMM:
		scaleHist("/mnt/e/SABREsim/det/kin3mc/kin3mc_7Li3He4He6Li2186keV_4He2H_7500keV_theta178218_phi_-2.125_2.125_detPlothistos.root", "1par/h1par_RecMissMassExE", "hMMM_RecoilExE_SIMSCALED", "/mnt/e/SABREsim/det/kin3mc/kin3mc_7Li3He4He6Li2186keV_4He2H_7500keV_theta178218_phi_-2.125_2.125_histos_MMM_Scaled.root", 6228./1229380.)
	
*/
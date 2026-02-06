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
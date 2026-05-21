#include <vector>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLeaf.h"


void DalitzPlots_9B_a5Li_vs_p8Be(TString filename_p8Be, TString filename_a5Li){

	TFile *infile_p8Be = new TFile(filename_p8Be, "READ");
	TFile *infile_a5Li = new TFile(filename_a5Li, "READ");

	if(!infile_p8Be || infile_p8Be->IsZombie() || !infile_a5Li || infile_a5Li->IsZombie()){
		std::cout << "eror with one or more file" << std::endl;
		return;
	}

	TTree *p8Be_tree = (TTree*)infile_p8Be->Get("InvMass_Mult3");
	TTree *a5Li_tree = (TTree*)infile_a5Li->Get("InvMass_Mult3");

	double val_p8Be[6];
	double val_a5Li[6];
	std::vector<TString> branches = {"012", "021", "102", "120", "201", "210"};

	for(int i=0; i<6; i++){
		//p8Be_tree->SetBranchAddress(branches[i]+".IM2_int", &val_p8Be[i]);
		//a5Li_tree->SetBranchAddress(branches[i]+".IM2_int", &val_a5Li[i]);
		TLeaf *leaf_p8Be = p8Be_tree->GetLeaf(branches[i], "IM2_int");
		TLeaf *leaf_a5Li = a5Li_tree->GetLeaf(branches[i], "IM2_int");

		if(!leaf_p8Be || !leaf_a5Li){
			std::cout << "Error, couldn't find leaf IM2_int" << std::endl;
			return;
		}

		leaf_p8Be->SetAddress(&val_p8Be);
		leaf_a5Li->SetAddress(&val_a5Li);
	}

	TFile *outfile = new TFile("DalitzPlots_9B_a5Li_vs_p8Be.root", "RECREATE");
	TH2D *DalitzPlot = new TH2D("DalitzPlot", "Dalitz Plot 9B, a+5Li vs p+8Be;IM2_int (8Be) MeV^2/c^4;IM2_int (5Li) MeV^2/c^4", 1520/8, 55.55e6, 55.92e6, 800/8, 21.75e6, 21.94e6);
	//TH2D *DalitzPlot = new TH2D("DalitzPlot", "Dalitz Plot 9B, a+5Li vs p+8Be;IM2_int (8Be) MeV^2/c^4;IM2_int (5Li) MeV^2/c^4", 1520, 0, 55.92e6, 800, 0, 21.94e6);

	Long64_t nentries = p8Be_tree->GetEntries();
	for(Long64_t i=0; i<nentries; i++){
		p8Be_tree->GetEntry(i);
		a5Li_tree->GetEntry(i);

		for(int x_idx=0; x_idx<6; x_idx++){
			for(int y_idx=0; y_idx<6; y_idx++){
				DalitzPlot->Fill(val_p8Be[x_idx], val_a5Li[y_idx]);
			}
		}
	}

	outfile->cd();
	DalitzPlot->Write();

	DalitzPlot->SetDirectory(gROOT);

	TCanvas *c = new TCanvas("c1", "Dalitz Plot", 1000, 1000);
	DalitzPlot->Draw("COLZ");

	outfile->Close();
	infile_p8Be->Close();
	infile_a5Li->Close();
}
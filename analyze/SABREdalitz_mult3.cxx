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
#include "TTree.h"
#include "TSystem.h"


//void SABREsim_DalitzPlots_9B_a5Li_vs_p8Be(int recExkeV, TString simDecayString){
void SABREsim_DalitzPlots_9B_a5Li_vs_p8Be(int recExkeV){

	TString filename_p8Be = Form("/mnt/e/SABREsim/analyze/may21/det/b10ha_7.5MeV_9B_ex%dkeV_tree_mult3_SABREanalyzed8Be.root", recExkeV);
	TString filename_a5Li = Form("/mnt/e/SABREsim/analyze/may21/det/b10ha_7.5MeV_9B_ex%dkeV_tree_mult3_SABREanalyzed5Li.root", recExkeV);

	TFile *infile_p8Be = new TFile(filename_p8Be, "READ");
	TFile *infile_a5Li = new TFile(filename_a5Li, "READ");

	if(!infile_p8Be || infile_p8Be->IsZombie() || !infile_a5Li || infile_a5Li->IsZombie()){
		std::cout << "eror with one or more file" << std::endl;
		return;
	}

	TTree *p8Be_tree = (TTree*)infile_p8Be->Get("InvMass_Mult3");
	TTree *a5Li_tree = (TTree*)infile_a5Li->Get("InvMass_Mult3");

	// double val_p8Be[6];
	// double val_a5Li[6];
	std::vector<TString> branches = {"012", "021", "102", "120", "201", "210"};

	TLeaf *leaves_p8Be[6];
	TLeaf *leaves_a5Li[6];

	for(int i=0; i<6; i++){
		//p8Be_tree->SetBranchAddress(branches[i]+".IM2_int", &val_p8Be[i]);
		//a5Li_tree->SetBranchAddress(branches[i]+".IM2_int", &val_a5Li[i]);
		leaves_p8Be[i] = p8Be_tree->GetLeaf(branches[i], "IM2_int");
		leaves_a5Li[i] = a5Li_tree->GetLeaf(branches[i], "IM2_int");

		if(!leaves_p8Be[i] || !leaves_a5Li[i]){
			std::cout << "Error, couldn't find leaf IM2_int" << std::endl;
			return;
		}

		// leaves_p8Be[i]->SetAddress(&val_p8Be[i]);
		// leaves_a5Li->SetAddress(&val_a5Li[i]);
	}

	TString dirPath = gSystem->DirName(filename_p8Be);
	TString baseOutfilename = Form("DalitzPlots_9Bex%dkeV.root", recExkeV);
	TString outfilename = dirPath+"/"+baseOutfilename;
	TFile *outfile = new TFile(outfilename, "RECREATE");

	TString xtitle = Form("IM2_int (8Be) MeV^2/c^4");
	TString ytitle = Form("IM2_int (5Li) MeV^2/c^4");
	//TString legibleDecayString = simDecayString.EqualTo("p8Be") ? "p + 8Be (gs)" : "a + 5Li (gs)" ;
	TString htitle = Form("Dalitz Plot 9B (%d keV);%s;%s", recExkeV, xtitle.Data(), ytitle.Data());

	TH2D *DalitzPlot = new TH2D("DalitzPlot", htitle, 1520/4, 55.55e6, 55.92e6, 800/4, 21.75e6, 21.94e6);
	//TH2D *DalitzPlot = new TH2D("DalitzPlot", "Dalitz Plot 9B, a+5Li vs p+8Be;IM2_int (8Be) MeV^2/c^4;IM2_int (5Li) MeV^2/c^4", 1520, 0, 55.92e6, 800, 0, 21.94e6);

	Long64_t nentries = p8Be_tree->GetEntries();
	for(Long64_t i=0; i<nentries; i++){
		p8Be_tree->GetEntry(i);
		a5Li_tree->GetEntry(i);

		for(int x_idx=0; x_idx<6; x_idx++){
			for(int y_idx=0; y_idx<6; y_idx++){
				DalitzPlot->Fill(leaves_p8Be[x_idx]->GetValue(), leaves_a5Li[y_idx]->GetValue());
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

//only difference between below and above is naming scheme of output and not paassing recExKeV to the data one
void exp_DalitzPlots_9B_a5Li_vs_p8Be(TString filename_p8Be, TString filename_a5Li){

	TFile *infile_p8Be = new TFile(filename_p8Be, "READ");
	TFile *infile_a5Li = new TFile(filename_a5Li, "READ");

	if(!infile_p8Be || infile_p8Be->IsZombie() || !infile_a5Li || infile_a5Li->IsZombie()){
		std::cout << "eror with one or more file" << std::endl;
		return;
	}

	TTree *p8Be_tree = (TTree*)infile_p8Be->Get("InvMass_Mult3");
	TTree *a5Li_tree = (TTree*)infile_a5Li->Get("InvMass_Mult3");

	// double val_p8Be[6];
	// double val_a5Li[6];
	std::vector<TString> branches = {"012", "021", "102", "120", "201", "210"};

	TLeaf *leaves_p8Be[6];
	TLeaf *leaves_a5Li[6];

	for(int i=0; i<6; i++){
		//p8Be_tree->SetBranchAddress(branches[i]+".IM2_int", &val_p8Be[i]);
		//a5Li_tree->SetBranchAddress(branches[i]+".IM2_int", &val_a5Li[i]);
		leaves_p8Be[i] = p8Be_tree->GetLeaf(branches[i], "IM2_int");
		leaves_a5Li[i] = a5Li_tree->GetLeaf(branches[i], "IM2_int");

		if(!leaves_p8Be[i] || !leaves_a5Li[i]){
			std::cout << "Error, couldn't find leaf IM2_int" << std::endl;
			return;
		}

		// leaves_p8Be[i]->SetAddress(&val_p8Be);
		// leaf_a5Li->SetAddress(&val_a5Li);
	}

	TFile *outfile = new TFile("DalitzPlots_B10ha_3par_exp.root", "RECREATE");
	TH2D *DalitzPlot = new TH2D("DalitzPlot", "Dalitz Plot 9B, a+5Li vs p+8Be;IM2_int (8Be) MeV^2/c^4;IM2_int (5Li) MeV^2/c^4", 1520/8, 55.55e6, 55.92e6, 800/8, 21.75e6, 21.94e6);
	//TH2D *DalitzPlot = new TH2D("DalitzPlot", "Dalitz Plot 9B, a+5Li vs p+8Be;IM2_int (8Be) MeV^2/c^4;IM2_int (5Li) MeV^2/c^4", 1520, 0, 55.92e6, 800, 0, 21.94e6);

	Long64_t nentries = p8Be_tree->GetEntries();
	for(Long64_t i=0; i<nentries; i++){
		p8Be_tree->GetEntry(i);
		a5Li_tree->GetEntry(i);

		for(int x_idx=0; x_idx<6; x_idx++){
			for(int y_idx=0; y_idx<6; y_idx++){
				DalitzPlot->Fill(leaves_p8Be[x_idx]->GetValue(), leaves_a5Li[y_idx]->GetValue());
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
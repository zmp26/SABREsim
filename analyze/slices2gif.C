#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TList.h"
#include "TStyle.h"
#include "TString.h"



//BELOW PLOTS DALITZ BOUNDARIES ON TOP OF DALITZ PLOTS
void slices2gif(const char* filename, const char* dirname){
	//filename is path to ROOT file containing the Dalitz/Catania slices output from SABREanalyze_mult3.cxx::B10ha_SABREPID()/_Mult2()
	//dirname is ROOT directory for histograms of interest

	TString gifFilename = "slices_animation.gif";
	// TString rmgif = Form("rm %s", gifFilename.Data());
	// gSystem->Exec(rmgif.Data());
	gSystem->Unlink(gifFilename.Data());


	TFile *file = TFile::Open(filename, "READ");
	if(!file || file->IsZombie()){
		std::cerr << "Error: Could not open file " << filename << std::endl;
		return;
	}

	TDirectory *dir = file->GetDirectory(dirname);
	if(!dir){
		std::cerr << "Error: Directory " << dirname << " not found in file.";
		file->Close();
		return;
	}

	gSystem->mkdir("./png",kTRUE);

	TList *keylist = dir->GetListOfKeys();
	if(!keylist || keylist->GetSize() == 0){
		std::cerr << "No keys found in directory: " << dirname << std::endl;
		file->Close();
		delete file;
		return;
	}

	//lambda to extract index i from names such as:
	//hDalitz_17_1.7_1.8 --> i = 17
	//gDalitzBoundary_17 --> i = 17
	auto GetSliceIndex = [](const TString& name) -> int {
		Ssiz_t firstunderscore = name.First("_");

		if(firstunderscore == kNPOS){
			return -1;
		}

		TString suffix = name;
		suffix.Remove(0, firstunderscore+1);

		//Atoi() stops at the next underscore/decimal point/etc
		return suffix.Atoi();
	};

	std::vector<TH2D*> histograms;

	TIter next(keylist);
	TKey *key = nullptr;
	while((key = (TKey*)next())){
		TObject *obj = key->ReadObj();
		if(obj && obj->InheritsFrom(TH2D::Class())){
			TH2D *h = static_cast<TH2D*>(obj);
			TString hname = h->GetName();
			if(hname.BeginsWith("hDalitz_")) histograms.push_back(h);
		}
	}

	std::sort(histograms.begin(), histograms.end(), [&](TH2D* a, TH2D* b) {
		return GetSliceIndex(a->GetName()) < GetSliceIndex(b->GetName());
	});

	if(histograms.empty()){
		std::cerr << "No 2D histograms found in " << dirname << std::endl;
		file->Close();
		return;
	}

	TCanvas *canvas = new TCanvas("c_gif", "Slice Viewer", 800, 800);
	gStyle->SetOptStat(0);

	for(size_t i=0; i<histograms.size(); i++){
		TH2D *h = histograms[i];

		const int sliceIndex = GetSliceIndex(h->GetName());
		TString boundaryname = Form("gDalitzBoundary_%02d",sliceIndex);
		TGraph* boundary = dynamic_cast<TGraph*>(dir->Get(boundaryname));
		//TGraph* boundary = nullptr;

		canvas->cd();
		canvas->Clear();

		//h->SetMinimum(0);
		h->Draw("COLZ");

		if(boundary){
			boundary->SetLineColor(kRed+1);
			boundary->SetLineWidth(2);
			boundary->SetLineStyle(1);

			boundary->Draw("L SAME");//SAME draws on same canvas, L adds
		} else {
			std::cerr << "Warning: no boundary graph named " << boundaryname << " found for histogram " << h->GetName() << std::endl;
		}

		canvas->Modified();
		canvas->Update();

		TString pngPath = Form("./png/%s.png", h->GetName());
		canvas->SaveAs(pngPath);

		if(i == histograms.size() - 1){
			canvas->Print(gifFilename + "++");
		} else {
			canvas->Print(gifFilename + "+50");
		}
	}

	std::cout << "\nFinished!" << std::endl;
	std::cout << "Saved " << histograms.size() << " PNG file to ./png/" << std::endl;
	std::cout << "Created animated gif: " << gifFilename << std::endl;

	delete canvas;
	file->Close();
	delete file;

}



/*
//below plots dalitz/catania by themselves

void slices2gif(const char* filename, const char* dirname){
	//filename is path to ROOT file containing the Dalitz/Catania slices output from SABREanalyze_mult3.cxx::B10ha_SABREPID()/_Mult2()
	//dirname is ROOT directory for histograms of interest

	TString gifFilename = "slices_animation.gif";
	TString rmgif = Form("rm %s", gifFilename.Data());
	gSystem->Exec(rmgif.Data());


	TFile *file = TFile::Open(filename, "READ");
	if(!file || file->IsZombie()){
		std::cerr << "Error: Could not open file " << filename << std::endl;
		return;
	}

	TDirectory *dir = file->GetDirectory(dirname);
	if(!dir){
		std::cerr << "Error: Directory " << dirname << " not found in file.";
		file->Close();
		return;
	}

	gSystem->mkdir("./png",kTRUE);

	TList *keylist = dir->GetListOfKeys();
	if(!keylist || keylist->GetSize() == 0){
		std::cerr << "No keys found in directory: " << dirname << std::endl;
		file->Close();
		return;
	}

	std::vector<TH2D*> histograms;

	TIter next(keylist);
	TKey *key = nullptr;
	while((key = (TKey*)next())){
		TObject *obj = key->ReadObj();
		if(obj && obj->InheritsFrom(TH2::Class())){
			histograms.push_back((TH2D*)obj);
		}
	}

	std::sort(histograms.begin(), histograms.end(), [](TH2D* a, TH2D* b) {
		TString nameA = a->GetName();
		TString nameB = b->GetName();

		int indexA = TString(nameA.Data() + nameA.Index("_") + 1).Atoi();
		int indexB = TString(nameB.Data() + nameB.Index("_") + 1).Atoi();

		return indexA < indexB;
	});

	if(histograms.empty()){
		std::cerr << "No 2D histograms found in " << dirname << std::endl;
		file->Close();
		return;
	}

	TCanvas *canvas = new TCanvas("c_gif", "Slice Viewer", 800, 800);
	gStyle->SetOptStat(0);

	for(size_t i=0; i<histograms.size(); i++){
		TH2D *h = histograms[i];

		canvas->Clear();
		h->Draw("COLZ");

		TString pngPath = Form("./png/%s.png", h->GetName());
		canvas->SaveAs(pngPath);

		if(i == histograms.size() - 1){
			canvas->Print(gifFilename + "+20");
		} else {
			canvas->Print(gifFilename + "+50");
		}
	}

	std::cout << "\nFinished!" << std::endl;
	std::cout << "Saved " << histograms.size() << " PNG file to ./png/" << std::endl;
	std::cout << "Created animated gif: " << gifFilename << std::endl;

	delete canvas;
	file->Close();
	delete file;

}
*/
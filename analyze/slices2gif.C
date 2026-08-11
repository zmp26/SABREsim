#include <iostream>
#include <vector>
#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TList.h"

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
			canvas->Print(gifFilename + "+500");
		} else {
			canvas->Print(gifFilename + "+200");
		}
	}

	std::cout << "\nFinished!" << std::endl;
	std::cout << "Saved " << histograms.size() << " PNG file to ./png/" << std::endl;
	std::cout << "Created animated gif: " << gifFilename << std::endl;

	delete canvas;
	file->Close();
	delete file;

}
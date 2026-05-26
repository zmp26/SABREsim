#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"

void EffCor_9B_a5Li(){
	if(!gFile){
		std::cerr << "Error: No open ROOT file in this session!" << std::endl;
		return;
	}

	TH1D* hist = (TH1D*)(gFile->Get("px"));
	if(!hist){
		std::cerr << "Histo not found!" << std::endl;
		return;
	}

	TString txtFileName = "../efficiencies/9B_4He_p_eff_1690_6000_SABREres50keV_thresh100keV.txt";
	std::vector<std::pair<double,double>> threeParEffs;
	std::ifstream efffile(txtFileName.Data());
	if(!efffile.is_open()){
		std::cerr << "Error: failed to open efficiency file" << std::endl;
		return;
	}

	std::string line;
	while(std::getline(efffile, line)){
		if(line.empty() || line[0] == 'E') continue;

		std::stringstream ss(line);
		double recEx, bu1, bu2, bu3, bu12, bu23, bu13, bu123;
		while(ss >> recEx >> bu1 >> bu2 >> bu3 >> bu12 >> bu23 >> bu13 >> bu123){
			threeParEffs.push_back({recEx, bu123/100.});
		}
	}

	efffile.close();
	TGraph* effGraph = new TGraph();

	for(int i=0; i<(int)threeParEffs.size(); i++){
		effGraph->SetPoint(i, threeParEffs.at(i).first, threeParEffs.at(i).second);
	}

	if(effGraph->GetN() == 0){
		std::cerr << "Error: failed to read data from " << txtFileName << std::endl;
		delete effGraph;
		return;
	}

	int nPoints = effGraph->GetN();
	for(int i=0; i<nPoints; i++){
		double energykeV = effGraph->GetPointX(i);
		effGraph->SetPointX(i, energykeV/1000.);
	}

	std::cout << "Loaded " << nPoints << " efficiency points!" << std::endl;
	std::cout << "Efficiency correcting..." << std::endl;

	int nbins = hist->GetNbinsX();
	int correctedCount = 0;

	for(int bin=1; bin<=nbins; bin++){

		double binCenterMeV = hist->GetBinCenter(bin);

		double efficiency = effGraph->Eval(binCenterMeV);

		if(efficiency <= 0.){
			continue;
		}

		double content = hist->GetBinContent(bin);
		double error = hist->GetBinError(bin);

		hist->SetBinContent(bin, content / efficiency);
		hist->SetBinError(bin, error / efficiency);

		correctedCount++;

	}

	std::cout << "Done! corrected " << correctedCount << " entries" << std::endl;

	//hist->Draw("HIST");

	TString oldFileName = gFile->GetName();
	TString newFileName = oldFileName;

	if(newFileName.EndsWith(".root")){
		newFileName.ReplaceAll(".root", "_effCor.root");
	} else {
		newFileName.Append("_effCor.root");
	}

	TFile *outfile = new TFile(newFileName, "RECREATE");
	if(outfile && !outfile->IsZombie()){
		hist->Write();
		outfile->Close();
		delete outfile;
	} else {
		std::cerr << "Error: Could not create output file!" << std::endl;
 	}


	delete effGraph;
}
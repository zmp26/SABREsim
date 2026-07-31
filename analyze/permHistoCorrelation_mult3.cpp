#include "permHistoCorrelation_mult3.h"
#include <iostream>

permHistoCorrelation_mult3::permHistoCorrelation_mult3(const std::vector<TString>& permNames, TDirectory* targetDir) 
	: fPermNames(permNames) 
{
	// Create or enter sub-directory for neatness in ROOT output file
	TDirectory* corrDir = targetDir->mkdir("Chi2_CornerPlots");
	if(corrDir) corrDir->cd();
	else targetDir->cd();

	int nPerms = fPermNames.size();
	int nBins = 60;
	double minChi2 = 0.0;
	double maxChi2 = 15.0; // Adjust max chi2 display limit as needed

	for(int i = 0; i < nPerms; ++i){
		// 1D Marginal on diagonal
		TString h1Key = Form("chi2_%s", fPermNames[i].Data());
		TString h1Title = Form("#chi^{2} Distribution for %s; #chi^{2}_{%s}; Counts", 
								fPermNames[i].Data(), fPermNames[i].Data());
		h1Map[h1Key] = new TH1D(h1Key, h1Title, nBins, minChi2, maxChi2);

		// 2D Cross-correlations (Lower triangle: i > j)
		for(int j = 0; j < i; ++j){
			TString h2Key = Form("chi2_%s_vs_%s", fPermNames[i].Data(), fPermNames[j].Data());
			TString h2Title = Form("#chi^{2} %s vs %s; #chi^{2}_{%s}; #chi^{2}_{%s}", 
									fPermNames[i].Data(), fPermNames[j].Data(), 
									fPermNames[j].Data(), fPermNames[i].Data());

			h2Map[h2Key] = new TH2D(h2Key, h2Title, nBins, minChi2, maxChi2, nBins, minChi2, maxChi2);
		}

		//(1D permutation - 012) (permutation minus correct permutation --> greater than 0 means worse pick than 012, lesser than 0 means better pick than 012)
		TString difkey = Form("chi2_dif_%s_minus_012", fPermNames[i].Data());
		TString diftitle = Form("%s - 012;#chi^{2}_{%s} - #chi^{2}_{012}", fPermNames[i].Data(), fPermNames[i].Data());
		difMap[difkey] = new TH1D(difkey, diftitle, nBins*8, -maxChi2, maxChi2);
	}
}

void permHistoCorrelation_mult3::FillCorner(const std::vector<double>& chi2Values){
	if(chi2Values.size() != fPermNames.size()){
		std::cerr << "Error: chi2Values vector size mismatch in FillCorner!\n";
		return;
	}

	int nPerms = fPermNames.size();

	for(int i = 0; i < nPerms; ++i){
		// Fill 1D diagonal
		TString h1Key = Form("chi2_%s", fPermNames[i].Data());
		if(h1Map.count(h1Key)){
			h1Map[h1Key]->Fill(chi2Values[i]);
		}

		// Fill 2D off-diagonal (j < i)
		for(int j = 0; j < i; ++j){
			TString h2Key = Form("chi2_%s_vs_%s", fPermNames[i].Data(), fPermNames[j].Data());
			if(h2Map.count(h2Key)){
				// X-axis: Permutation j, Y-axis: Permutation i
				h2Map[h2Key]->Fill(chi2Values[j], chi2Values[i]);
			}
		}

		//Fill dif:
		TString difkey = Form("chi2_dif_%s_minus_012", fPermNames[i].Data());
		if(difMap.count(difkey)){
			difMap[difkey]->Fill(chi2Values[i] - chi2Values[0]);
		}
	}
}
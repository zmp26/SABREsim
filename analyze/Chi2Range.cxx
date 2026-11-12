#include <TH1.h>
#include <TMath.h>
#include <iostream>

//takes files and histo name in, good for calling on existing files
double Chi2Range(TString file1, TString histo1, TString file2, TString histo2, double xmin, double xmax, int& ndf, bool scale=true){
	TFile *infile1 = TFile::Open(file1, "READ");
	TFile *infile2 = TFile::Open(file2, "READ");
	if(!infile1 || !infile2){
		std::cout << "can't open at least one of the files" << std::endl;
		ndf = 0;
		return -1.;
	}

	TH1D *h1 = (TH1D*)infile1->Get(histo1);
	TH1D *h2 = (TH1D*)infile2->Get(histo2);
	if (!h1 || !h2) {
		std::cerr << "Error: One or both histogram pointers are null!" << std::endl;
		ndf = 0;
		return -1.;
	}

	// Find the corresponding bin indices for the given x range
	int binMin = h1->FindBin(xmin);
	int binMax = h1->FindBin(xmax);

	// Ensure bins are within valid histogram limits (1 to NbinsX)
	binMin = std::max(1, binMin);
	binMax = std::min(h1->GetNbinsX(), binMax);

	//get integral over specified range:
	double integral1 = h1->Integral(xmin, xmax);
	double integral2 = h2->Integral(xmin, xmax);

	TH1D *h1_scaled = (TH1D*)h1->Clone("h1_scaled");
	TH1D *h2_scaled = (TH1D*)h2->Clone("h2_scaled");

	if(scale && integral1 > 0 && integral2 > 0){
		if(integral1 > integral2){
			double scalefactor = integral2 / integral1;
			h1_scaled->Scale(scalefactor);
			std::cout << "Scaling down histogram 1 (" << histo1 << ") by factor " << scalefactor << std::endl;
		} else if(integral2 < integral1){
			double scalefactor = integral1 / integral2;
			h2_scaled->Scale(scalefactor);
			std::cout << "Scaling down histogram 2 (" << histo2 << ") by factor " << scalefactor << std::endl;
		}
	}

	double chi2 = 0.0;
	ndf = 0;

	for (int i = binMin; i <= binMax; ++i) {
		double y1 = h1_scaled->GetBinContent(i);
		double y2 = h2_scaled->GetBinContent(i);
		double e1 = h1_scaled->GetBinError(i);
		double e2 = h2_scaled->GetBinError(i);

		// Total variance for independent measurements: sigma_tot^2 = sigma_1^2 + sigma_2^2
		double errorSq = (e1 * e1) + (e2 * e2);

		// Skip bins with zero error to avoid division by zero
		if (errorSq > 0.0) {
			double diff = y1 - y2;
			chi2 += (diff * diff) / errorSq;
			ndf++;
		}
	}

	//clean up temp clone objects:
	delete h1_scaled;
	delete h2_scaled;

	return chi2;
}

//takes pointers in, good for using in-session
double Chi2Range(TH1* h1, TH1* h2, double xmin, double xmax, int& ndf, bool scale = true) {
	// Check if the input histogram pointers are valid
	if (!h1 || !h2) {
		std::cerr << "Error: One or both input histogram pointers are null!" << std::endl;
		ndf = 0;
		return -1.;
	}

	// Find the corresponding bin indices for the given x range
	int binMin = h1->FindBin(xmin);
	int binMax = h1->FindBin(xmax);

	// Ensure bins are within valid histogram limits (1 to NbinsX)
	binMin = std::max(1, binMin);
	binMax = std::min(h1->GetNbinsX(), binMax);

	// Get integral over specified range
	double integral1 = h1->Integral(xmin, xmax);
	double integral2 = h2->Integral(xmin, xmax);

	// Clone histograms so we don't modify the originals in memory
	TH1D *h1_scaled = (TH1D*)h1->Clone("h1_scaled");
	TH1D *h2_scaled = (TH1D*)h2->Clone("h2_scaled");

	if (scale && integral1 > 0 && integral2 > 0) {
		if (integral1 > integral2) {
			double scalefactor = integral2 / integral1;
			h1_scaled->Scale(scalefactor);
			std::cout << "Scaling down histogram 1 (" << h1->GetName() << ") by factor " << scalefactor << std::endl;
		} else if (integral2 > integral1) {
			double scalefactor = integral1 / integral2;
			h2_scaled->Scale(scalefactor);
			std::cout << "Scaling down histogram 2 (" << h2->GetName() << ") by factor " << scalefactor << std::endl;
		}
	}

	double chi2 = 0.0;
	ndf = 0;

	for (int i = binMin; i <= binMax; ++i) {
		double y1 = h1_scaled->GetBinContent(i);
		double y2 = h2_scaled->GetBinContent(i);
		double e1 = h1_scaled->GetBinError(i);
		double e2 = h2_scaled->GetBinError(i);

		// Total variance for independent measurements: sigma_tot^2 = sigma_1^2 + sigma_2^2
		double errorSq = (e1 * e1) + (e2 * e2);

		// Skip bins with zero error to avoid division by zero
		if (errorSq > 0.0) {
			double diff = y1 - y2;
			chi2 += (diff * diff) / errorSq;
			ndf++;
		}
	}

	// Clean up temporary clone objects to prevent memory leaks
	delete h1_scaled;
	delete h2_scaled;

	return chi2;
}
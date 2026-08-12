#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TDirectory.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TColor.h"
#include "TStyle.h"
#include "TString.h"

void makeHelicityDalitzSlices(const char* outFilename = "helicity_dalitz_slices.root"){
	TMassTable fMassTable;
	fMassTable.Init("/home/zachpurcell/masstable/masstable.dat");

	double mass_p = fMassTable.GetNuclearMassMeV("H", 1);
	double mass_a = fMassTable.GetNuclearMassMeV("He", 4);
	double massgs_9B = fMassTable.GetNuclearMassMeV("B", 9);

	auto Kallen = [](double a, double b, double c) { return a*a + b*b + c*c - 2.0*(a*b + a*c + b*c); };

	auto MakeDalitzBoundary = [&](double W, int sliceIndex, double sliceEMin, double sliceEMax, int nPoints = 500) -> TGraph*
	{
		double xmin = 4.0*mass_a*mass_a;
		double xmax = std::pow(W - mass_p, 2);

		std::vector<double> x;
		std::vector<double> y;

		x.reserve(2*nPoints);
		y.reserve(2*nPoints);

		for(int j = 0; j < nPoints; ++j) {
			double frac = static_cast<double>(j)/(nPoints - 1);

			double sAA = xmin + frac*(xmax - xmin);

			double lambdaAA = Kallen(sAA, mass_a*mass_a, mass_a*mass_a);

			double lambdaW = Kallen(W*W, sAA, mass_p*mass_p);

			double radical = std::sqrt(std::max(0.0, lambdaAA*lambdaW));

			double yCenter = mass_p*mass_p + mass_a*mass_a + 0.5*(W*W - mass_p*mass_p - sAA);

			double yPlus = yCenter + radical/(2.0*sAA);

			x.push_back(sAA);
			y.push_back(yPlus);
		}

		for(int j = nPoints - 1; j >= 0; j--) {
			double frac = static_cast<double>(j)/(nPoints - 1);

			double sAA = xmin + frac*(xmax - xmin);

			double lambdaAA = Kallen(sAA, mass_a*mass_a, mass_a*mass_a);

			double lambdaW = Kallen(W*W, sAA, mass_p*mass_p);

			double radical = std::sqrt(std::max(0.0, lambdaAA*lambdaW));

			double yCenter = mass_p*mass_p + mass_a*mass_a + 0.5*(W*W - mass_p*mass_p - sAA);

			double yMinus = yCenter - radical/(2.0*sAA);

			x.push_back(sAA);
			y.push_back(yMinus);
		}

		TGraph* boundary = new TGraph(
			static_cast<int>(x.size()), x.data(), y.data()
		);

		boundary->SetName(
			Form("gDalitzBoundary_%02d", sliceIndex)
		);

		boundary->SetTitle(Form("Dalitz boundary, E_{x}=%.1f-%.1f MeV;M^{2}_{#alpha#alpha};M^{2}_{p#alpha}",sliceEMin, sliceEMax));

		boundary->SetLineColor(kBlack);
		boundary->SetLineWidth(2);
		boundary->SetLineStyle(1);

		return boundary;
	};

	double eMin = 0.0;
	double eMax = 7.0;
	double eStep = 0.1;

	int nSlices = static_cast<int>((eMax - eMin)/eStep);

	TFile* outfile = new TFile(outFilename, "RECREATE");

	if(!outfile || outfile->IsZombie()) {
		std::cerr << "Could not create " << outFilename << std::endl;
		return;
	}

	TDirectory* sliceDir = outfile->mkdir("slices");

	//histogram names intentionally match hDalitz_* so that slices2gif macro can read them.
	for(int i = 0; i < nSlices; ++i) {
		double currentEMin = eMin + i*eStep;
		double currentEMax = eMin + (i + 1)*eStep;
		double currentECenter = 0.5*(currentEMin + currentEMax);

		//assumes Ex is Ex(9B)
		double W = massgs_9B + currentECenter;

		TString name = Form("hDalitz_%02d_%.1f_%.1f", i, currentEMin, currentEMax);

		TString title = Form("cos(#theta^{h}_{#alpha}) over allowed Dalitz space, E_{x}=%.1f-%.1f MeV;M^{2}_{#alpha#alpha} [MeV^{2}];M^{2}_{p#alpha} [MeV^{2}]", currentEMin, currentEMax);

		TH2D* hCosTheta = new TH2D(name.Data(), title.Data(), 240, 55.57e6, 55.69e6, 160, 21.76e6, 21.84e6);

		hCosTheta->SetDirectory(sliceDir);

		hCosTheta->SetMinimum(-1.0);
		hCosTheta->SetMaximum(+1.0);

		for(int xBin = 1; xBin <= hCosTheta->GetNbinsX(); ++xBin) {
			double sAA = hCosTheta->GetXaxis()->GetBinCenter(xBin);

			double lambdaAA = Kallen(sAA, mass_a*mass_a, mass_a*mass_a);

			double lambdaW = Kallen(W*W, sAA, mass_p*mass_p);

			if(lambdaAA <= 0.0 || lambdaW <= 0.0) continue;

			double deltaY = std::sqrt(lambdaAA*lambdaW)/(2.0*sAA);

			if(deltaY <= 1.0e-12) {
				continue;
			}

			double yCenter =
				mass_p*mass_p + mass_a*mass_a +
				0.5*(W*W - mass_p*mass_p - sAA);

			double yMinus = yCenter - deltaY;
			double yPlus  = yCenter + deltaY;

			for(int yBin = 1; yBin <= hCosTheta->GetNbinsY(); yBin++) {
				double sPA = hCosTheta->GetYaxis()->GetBinCenter(yBin);

				if(sPA < yMinus || sPA > yPlus) {
					continue;
				}

				double cosThetaH = (yCenter - sPA)/deltaY;
				cosThetaH = std::max(-1.0, std::min(1.0, cosThetaH));
				hCosTheta->SetBinContent(xBin, yBin, cosThetaH);
			}
		}

		TGraph* boundary = MakeDalitzBoundary(W, i, currentEMin, currentEMax);

		sliceDir->cd();
		boundary->Write();
	}

	outfile->cd();
	outfile->Write();
	outfile->Close();

	delete outfile;

	std::cout << "Created helicity-coordinate Dalitz slices in "
			  << outFilename << std::endl;
}
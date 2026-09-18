#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TH2D.h>
#include <TPaletteAxis.h>
#include <iostream>

void drawDalitzBoundaries(const char* filename) {
	TFile *file = TFile::Open(filename, "READ");
	if (!file || file->IsZombie()) {
		std::cerr << "Error opening file!" << std::endl;
		return;
	}

	TCanvas *canvas = new TCanvas("c1", "Dalitz Boundaries", 600, 600);
	
	// Set right margin to leave room for the color bar
	canvas->SetRightMargin(0.15);

	gStyle->SetPalette(kRainbow);
	Int_t nColors = gStyle->GetNumberOfColors();
	const Int_t totalGraphs = 350;

	// Define physical energy limits (adjust these to your actual energy values)
	Double_t minEnergy = 0.0;  // GeV or MeV
	Double_t maxEnergy = 7.0;  // GeV or MeV

	// 1. Draw a dummy 2D histogram off-screen or underneath to generate the palette axis
	TH2D *hDummy = new TH2D("hDummy", "Dalitz Boundaries", 100, 55.57e6, 55.69e6, 100, 21.76e6, 21.84e6);
	hDummy->SetMinimum(minEnergy);
	hDummy->SetMaximum(maxEnergy);
	hDummy->SetStats(0); // Disable stats box
	hDummy->GetZaxis()->SetTitle("E_{CM}(^{9}B) [MeV]"); // Label for the color bar
	hDummy->Draw("COLZ"); // "COLZ" forces creation of TPaletteAxis

	// 2. Draw all TGraph objects on top using palette color index
	for (Int_t i = 0; i < totalGraphs; ++i) {
		TString graphPath = TString::Format("slices/Dalitz/gDalitzBoundary_%03d", i);
		TGraph *graph = (TGraph*)file->Get(graphPath);

		// if (!graph){
		// 	std::cout << graphPath << " not found, continuing..." << std::endl;
		// 	continue;
		// } else {
		// 	std::cout << graphPath << "found!" << std::endl;
		// }

		// Calculate mapped color index
		Int_t colorIdx = gStyle->GetColorPalette((i * (nColors - 1)) / (totalGraphs - 1));

		graph->SetLineColor(colorIdx);
		graph->SetLineWidth(2);
		
		// Overlay onto the dummy histogram frame
		graph->Draw("L SAME");
	}

	// 3. Redraw the canvas and adjust the color bar position if needed
	canvas->Update();

	// Access the generated TPaletteAxis object to tweak palette formatting
	TPaletteAxis *palette = (TPaletteAxis*)hDummy->GetListOfFunctions()->FindObject("palette");
	if (palette) {
		palette->SetX1NDC(0.86);
		palette->SetX2NDC(0.89);
		palette->SetY1NDC(0.10);
		palette->SetY2NDC(0.90);
		canvas->Modified();
	}

	canvas->Draw();
}
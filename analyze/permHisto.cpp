#include "permHisto.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <vector>

permHisto::permHisto(TString permName, TDirectory* targetDir){
	targetDir->cd();

	std::vector<TString> particles = {"intermediate", "frag1", "frag2", "frag3"};

	//invariant mass and ex energy histograms
	Register1D(permName, "intermediateIM", "Intermediate Invariant Mass;MeV/c^{2}", 29000, 4600, 7500);
	Register1D(permName, "intermediateEx", "Intermediate Ex;MeV", 525, -1, 20);
	Register1D(permName, "RecoilEx", "Recoil Ex;MeV", 525, -1, 20);

	for(const auto& p : particles){
		Register1D(permName, p+"vcm_meas", p+" Velocity CM (meas);c", 250, 0, 0.25);
		Register1D(permName, p+"vcm_expect", p+" Velocity CM (expect);c", 250, 0, 0.25);
		Register1D(permName, p+"vcm_delta", p+" Velocity CM (meas - expect);c", 500, -0.25, 0.25);

		Register1D(permName, p+"kecm_meas", p+" KE CM (meas); MeV", 500, 0, 5);
		Register1D(permName, p+"kecm_expect", p+" KE CM (expect);MeV", 500, 0, 5);
		Register1D(permName, p+"kecm_delta", p+" KE CM (meas - expect);MeV", 1000, -5, 5);

		Register1D(permName, p+"thetacm", p+" Theta CM; deg", 36, 0, 180);
		Register1D(permName, p+"phicm", p+" Phi CM; deg", 72, 0, 360);
		Register2D(permName, p+"thetacmvsphicm", p+" ThetaCM vs PhiCM;deg;deg", 72, 0, 360, 36, 0, 180);
	}

	//decay1
	Register1D(permName, "ecm1_meas", "E_{cm} Decay 1 (meas);MeV", 600, -1, 5);
	Register1D(permName, "ecm1_expect", "E_{cm} Decay 1 (expect);MeV", 600, -1, 5);
	Register1D(permName, "ecm1_delta", "E_{cm} Decay 1 (meas - expect);MeV", 1000, -5, 5);
	Register2D(permName, "intermediatevcmVSfrag1vcm", "intermediate Vcm VS frag1 Vcm", 250, 0, 0.25, 250, 0, 0.25);
	Register2D(permName, "intermediatekecmVSfrag1kecm", "intermediate KEcm VS frag1 KEcm", 600, -1, 5, 600, -1, 5);

	//decay2
	Register1D(permName, "ecm2_meas", "E_{cm} Decay 2 (meas);MeV", 600, -1, 5);
	Register1D(permName, "ecm2_expect", "E_{cm} Decay 2 (expect);MeV", 600, -1, 5);
	Register1D(permName, "ecm2_delta", "E_{cm} Decay 2 (meas - expect);MeV", 1000, -5, 5);
	Register2D(permName, "frag2vcmVSfrag3vcm", "frag2 Vcm VS frag3 Vcm", 250, 0, 0.25, 250, 0, 0.25);
	Register2D(permName, "frag2kecmVSfrag3kecm", "frag2 KEcm VS frag3 KEcm", 600, -1, 5, 600, -1, 5);

	//sequential decay energies
	Register2D(permName, "ecm1VSecm2", "ECM1 vs ECM2; E_{CM} Decay 2 (MeV); E_{CM} Decay 1 (MeV)", 600, -1, 5, 600, -1 , 5);
	Register2D(permName, "ecm1deltaVSecm2delta", "ECM1 Delta vs ECM2 Delta; Decay 2 (MeV); Decay 1 (MeV)", 1000, -5, 5, 1000, -5, 5);

}

void permHisto::Register1D(TString name, TString key, TString title, int xbins, double xmin, double xmax){
	hMap[key] = new TH1D(name + "_" + key, title, xbins, xmin, xmax);
}

void permHisto::Register2D(TString name, TString key, TString title, int xbins, double xmin, double xmax, int ybins, int ymin, int ymax){
	hMap[key] = new TH2D(name + "_" + key, title, xbins, xmin, xmax, ybins, ymin, ymax);
}

void permHisto::Fill(TString key, double x){
	if(hMap.count(key)) hMap[key]->Fill(x);
	else std::cout << key << " not found in hMap!" << std::endl;
}

void permHisto::Fill(TString key, double x, double y){
	if(hMap.count(key)) hMap[key]->Fill(x, y);
	else std::cout << key << " not found in hMap!" << std::endl;
}
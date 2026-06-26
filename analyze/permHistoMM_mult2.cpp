#include "permHistoMM_mult2.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <vector>

permHistoMM_mult2::permHistoMM_mult2(TString permName, TDirectory* targetDir){
	targetDir->cd();

	std::vector<TString> particles = {"intermediate", "frag1", "frag2", "frag3"};

	//invariant mass and ex energy histograms
	//Register1D(permName, "intermediateIM", "Intermediate Invariant Mass;MeV/c^{2}", 5000, 0, 10000);
	Register1D(permName, "intermediateEx", "Intermediate Ex;MeV",300,-5,7);
	Register1D(permName, "RecoilEx", "Recoil Ex;MeV",300,-5,7);
	Register2D(permName, "RecoilEx_IMvsSPS", "Recoil Ex IM vs SPS;SPS MeV;IM MeV",200,-1,7,200,-1,7);
	Register1D(permName, "RecoilExDif", "Recoil Ex (SPS - IM);MeV", 100, -2, 2);
	Register2D(permName, "intermediateExIMvsSPS", "intermediate Ex IM vs SPS;SPS MeV;IM MeV",200,-1,7,200,-1,7);

	for(const auto& p : particles){
		Register1D(permName, p+"vcm_meas", p+" Velocity CM (meas);c", 5000, 0, 0.10);
		//Register1D(permName, p+"vcm_expect", p+" Velocity CM (expect);c", 5000, 0, 0.10);
		//Register1D(permName, p+"vcm_delta", p+" Velocity CM (meas - expect);c", 5000, -0.1, 0.1);
		Register2D(permName, p+"vcm_TransverseVSLongitudinal", p+" Velocity CM Transverse vs Longitudinal;c;c", 100, 0, 0.10, 100, 0, 0.10);

		Register1D(permName, p+"kecm_meas", p+" KE CM (meas); MeV", 500, 0, 5);
		//Register1D(permName, p+"kecm_expect", p+" KE CM (expect);MeV", 500, 0, 5);
		//Register1D(permName, p+"kecm_delta", p+" KE CM (meas - expect);MeV", 1000, -5, 5);

		Register1D(permName, p+"thetacm", p+" Theta CM; deg", 180, 0, 180);
		Register1D(permName, p+"phicm", p+" Phi CM; deg", 360, 0, 360);
		Register2D(permName, p+"thetacmvsphicm", p+" ThetaCM vs PhiCM;deg;deg", 360, 0, 360, 180, 0, 180);

		Register2D(permName, p+"vcmVSthetacm", p+" V CM vs Theta CM", 180, 0, 180, 100, 0, 0.1);
		Register2D(permName, p+"kecmVSthetacm", p+" KE CM vs Theta CM", 180, 0, 180, 500, 0, 5);
	}

	//missing mass specific histograms (note that these will change meaning between permutations!)
	Register1D(permName, "MissingMass", "Reconstructed Missing Mass;MeV/c^{2}", 200, 0, 4000);
	Register2D(permName, "intermediateExVSMissingMass", "intermediate vs ExMissing Mass;MeV;MeV/c^{2}", 500, 0, 10000, 200, -1, 7);

	//decay1
	Register1D(permName, "ecm1_meas", "E_{cm} Decay 1 (meas);MeV", 150, -1, 5);
	//Register1D(permName, "ecm1_expect", "E_{cm} Decay 1 (expect);MeV", 150, -1, 5);
	//Register1D(permName, "ecm1_delta", "E_{cm} Decay 1 (meas - expect);MeV", 250, -5, 5);
	Register2D(permName, "ecm1measVSintermediatethetacm", "E_{cm} Decay 1 (meas) vs Intermediate ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
	Register2D(permName, "ecm1measVSfrag1thetacm", "E_{cm} Decay 1 (meas) vs frag1 ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
	Register2D(permName, "ecm1measVSfrag2thetacm", "E_{cm} Decay 1 (meas) vs frag2 ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
	Register2D(permName, "ecm1measVSfrag3thetacm", "E_{cm} Decay 1 (meas) vs frag3 ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
	Register1D(permName, "decay1_VCM", "VCM Decay 1;c", 1000, 0, 0.1);
	Register2D(permName, "decay1_VCM_TransverseVSLongitudinal", "Decay 1 VCM Transverse Vs Longitudinal", 100, 0, 0.10, 100, 0, 0.10);
	Register1D(permName, "decay1_thetaCMsum", "decay1 Theta CM Sum", 80, 170, 190);
	Register1D(permName, "decay1_phiCMdiff", "decay1 Phi CM Diff", 80, 170, 190);
	Register1D(permName, "decay1_relLabAngle", "decay1 Relative Lab Angle", 360, 0, 180);
	Register2D(permName, "intermediatevcmVSfrag1vcm", "intermediate Vcm VS frag1 Vcm", 1000, 0, 0.1, 1000, 0, 0.1);
	Register2D(permName, "intermediatekecmVSfrag1kecm", "intermediate KEcm VS frag1 KEcm", 600, -1, 5, 600, -1, 5);

	//decay2
	Register1D(permName, "ecm2_meas", "E_{cm} Decay 2 (meas);MeV", 150, -1, 5);
	//Register1D(permName, "ecm2_expect", "E_{cm} Decay 2 (expect);MeV", 150, -1, 5);
	//Register1D(permName, "ecm2_delta", "E_{cm} Decay 2 (meas - expect);MeV", 250, -5, 5);
	Register2D(permName, "ecm2measVSintermediatethetacm", "E_{cm} Decay 2 (meas) vs Intermediate ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
	Register2D(permName, "ecm2measVSfrag1thetacm", "E_{cm} Decay 2 (meas) vs frag1 ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
	Register2D(permName, "ecm2measVSfrag2thetacm", "E_{cm} Decay 2 (meas) vs frag2 ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
	Register2D(permName, "ecm2measVSfrag3thetacm", "E_{cm} Decay 2 (meas) vs frag3 ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
	Register1D(permName, "decay2_VCM", "VCM Decay 2;c", 1000, 0, 0.1);
	Register2D(permName, "decay2_VCM_TransverseVSLongitudinal", "Decay 2 VCM Transverse Vs Longitudinal",  100, 0, 0.10, 100, 0, 0.10);
	Register1D(permName, "decay2_thetaCMsum", "decay2 Theta CM Sum", 80, 170, 190);
	Register1D(permName, "decay2_phiCMdiff", "decay2 Phi CM Diff", 80, 170, 190);
	Register1D(permName, "decay2_relLabAngle", "decay2 Relative Lab Angle", 360, 0, 180);
	Register2D(permName, "frag2vcmVSfrag3vcm", "frag2 Vcm VS frag3 Vcm", 1000, 0, 0.1, 1000, 0, 0.1);
	Register2D(permName, "frag2kecmVSfrag3kecm", "frag2 KEcm VS frag3 KEcm", 600, -1, 5, 600, -1, 5);

	Register2D(permName, "decay2VSdecay1_relLabAngle", "decay2 vs decay1 Relative Lab Angle", 360, 0, 180, 360, 0, 180);

	//sequential decay energies
	Register2D(permName, "ecm1VSecm2", "ECM1 vs ECM2; E_{CM} Decay 2 (MeV); E_{CM} Decay 1 (MeV)", 150, -1, 5, 150, -1 , 5);
	Register2D(permName, "ecm1deltaVSecm2delta", "ECM1 Delta vs ECM2 Delta; Decay 2 (MeV); Decay 1 (MeV)", 250, -5, 5, 250, -5, 5);

	//catania plot
	Register2D(permName, "catania_plot", "Catania Plot (p_{4}^{2}/(2m) vs E_{4}-Q)", 100, 0, 10, 100, 0, 10);

}

void permHistoMM_mult2::Register1D(TString name, TString key, TString title, int xbins, double xmin, double xmax){
	hMap[key] = new TH1D(name + "_" + key, title, xbins, xmin, xmax);
}

void permHistoMM_mult2::Register2D(TString name, TString key, TString title, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax){
	hMap[key] = new TH2D(name + "_" + key, title, xbins, xmin, xmax, ybins, ymin, ymax);
}

void permHistoMM_mult2::Fill(TString key, double x){
	if(hMap.count(key)) hMap[key]->Fill(x);
	else std::cout << key << " not found in hMap!" << std::endl;
}

void permHistoMM_mult2::Fill(TString key, double x, double y){
	if(hMap.count(key)) hMap[key]->Fill(x, y);
	else std::cout << key << " not found in hMap!" << std::endl;
}
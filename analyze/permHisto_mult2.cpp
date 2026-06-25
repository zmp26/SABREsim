#include "permHisto_mult2.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <vector>

permHisto_mult2::permHisto_mult2(TString permName, TDirectory* targetDir){
	targetDir->cd();

	std::vector<TString> particles = {"frag1", "frag2"};

	Register1D(permName, "RecoilIM", "Recoil Invariant Mass;MeV/c^{2}", 2900, 4600, 7500);
	Register1D(permName, "RecoilEx", "Recoil Ex;MeV", 300, -5, 7);
	Register2D(permName, "RecoilEx_IMvsSPS", "Recoil Ex IM vs SPS; SPS MeV; IM MeV", 200,-1,7,200,-1,7);
	Register1D(permName, "RecoilExDif", "Recoil Ex (SPS - IM); MeV", 100, -2, 2);

	for(const auto& p : particles){
		Register1D(permName, p+"vcm_meas", p+" Velocity CM (meas);c", 5000, 0, 0.1);
		//Register1D(permName, p+"vcm_expect", p+" Velocity CM (expect);c", 500, 0, 0.1);
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

	// Decay Kinematics
	Register1D(permName, "ecm_meas", "E_{cm} Decay (meas);MeV", 150, -1, 5);
	//Register1D(permName, "ecm_expect", "E_{cm} Decay (expect);MeV", 150, -1, 5);
	//Register1D(permName, "ecm_delta", "E_{cm} Decay (meas - expect);MeV", 250, -5, 5);
    
	Register2D(permName, "ecmmeasVSfrag1thetacm", "E_{cm} Decay (meas) vs frag1 ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
	Register2D(permName, "ecmmeasVSfrag2thetacm", "E_{cm} Decay (meas) vs frag2 ThetaCM;MeV;deg", 180, 0, 180, 150, -1, 5);
    
	Register1D(permName, "decay_VCM", "VCM Decay;c", 1000, 0, 0.1);
	Register2D(permName, "decay_VCM_TransverseVSLongitudinal", "Decay VCM Transverse Vs Longitudinal", 100, 0, 0.10, 100, 0, 0.10);
	Register1D(permName, "decay_thetaCMsum", "Decay Theta CM Sum", 80, 170, 190);
	Register1D(permName, "decay_phiCMdiff", "Decay Phi CM Diff", 80, 170, 190);
	Register1D(permName, "decay_relLabAngle", "Decay Relative Lab Angle", 360, 0, 180);
    
	Register2D(permName, "frag1vcmVSfrag2vcm", "frag1 Vcm VS frag2 Vcm", 1000, 0, 0.1, 1000, 0, 0.1);
	Register2D(permName, "frag1kecmVSfrag2kecm", "frag1 KEcm VS frag2 KEcm", 600, -1, 5, 600, -1, 5);
}

void permHisto_mult2::Register1D(TString name, TString key, TString title, int xbins, double xmin, double xmax){
	hMap[key] = new TH1D(name + "_" + key, title, xbins, xmin, xmax);
}

void permHisto_mult2::Register2D(TString name, TString key, TString title, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax){
	hMap[key] = new TH2D(name + "_" + key, title, xbins, xmin, xmax, ybins, ymin, ymax);
}

void permHisto_mult2::Fill(TString key, double x){
	if(hMap.count(key)) hMap[key]->Fill(x);
	else std::cout << key << " not found in hMap!" << std::endl;
}

void permHisto_mult2::Fill(TString key, double x, double y){
	if(hMap.count(key)) hMap[key]->Fill(x,y);
	else std::cout << key << " not found in hMap!" << std::endl;
}
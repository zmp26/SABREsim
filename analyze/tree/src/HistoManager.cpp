#include "HistoManager.h"

void HistoManager::Register1D(std::string key, std::string title, int bins, double min, double max){
	dir->cd();
	hMap[key] = std::make_unique<TH1D>(key.c_str(), title.c_str(), bins, min, max);
}

void HistoManager::Register2D(std::string key, std::string title, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax){
	dir->cd();
	hMap[key] = std::make_unique<TH2D>(key.c_str(), title.c_str(), xbins, xmin, xmax, ybins, ymin, ymax);
}

// void HistoManager::RegisterNodeHistograms(){
// 	std::string p = prefix + "_" + name;

// 	Register1D(p+"_thetacm", name+" Theta CM;deg", 360, 0, 180);
// 	Register1D(p+"_phicm", name+" Phi CM;deg", 720, 0, 720);
// 	Register1D(p+"_thetaCMVsPhiCM",name+" ThetaCM Vs PhiCM;deg;deg",360,0,180,720,0,360);

// 	Register1D(p+"_vcm_meas", name+" Vcm (meas);c", 5000, 0, 0.1);
// 	Register1D(p+"_kecm_meas", name+" KE CM (meas);MeV", 500, 0, 5);
// }

void HistoManager::RegisterHistograms(const KinematicNode* node, std::string perm){
	std::string n = node->name;
	std::string base = perm + "_" + n;

	// Register1D(perm, n+"_vcm", n+" Vcm (meas);c", 5000, 0, 0.1);
	// Register1D(perm, n+"_kecm", n+" KE CM; MeV", 500, 0, 5);
	// Register1D(perm, n+"_thetacm",n+" Theta CM; deg", 360, 0, 180);
	// Register2D(perm, n+"_thetaCMVsPhiCM",n+"ThetaCM vs PhiCM;deg;deg", 360, 0, 180, 720, 0, 360);
	Register1D(base+"_vcm", n+" Vcm (meas);c", 5000, 0, 0.1);
	Register1D(base+"_kecm", n+" KEcm;MeV",500,0,5);
	Register1D(base+"_thetacm", n+" ThetaCM;deg",360,0,180);
	Register1D(base+"_phicm", n+" PhiCM;deg",720,0,360);
	Register2D(base+"_thetaCMVsPhiCM", n+" ThetaCM vs PhiCM;deg;deg",720,0,360,360,0,180);


	if(!node->children.empty()){
		//Register1D(perm, n+"_ecm", n+" Ecm;MeV", 500, -1, 10);
		Register1D(base+"_ecm", n+" Ecm;MeV", 500, -1, 10);
		for(const auto& child : node->children){
			RegisterHistograms(child.get(), perm);
		}
	}

}

void HistoManager::Fill(std::string key, double x){
	if(hMap.count(key)) hMap[key]->Fill(x);
}

void HistoManager::Fill(std::string key, double x, double y){
	if(hMap.count(key)) hMap[key]->Fill(x,y);
}
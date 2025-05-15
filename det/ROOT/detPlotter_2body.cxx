/*

*/
using namespace std;
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TCutG.h"
#include "TRandom.h"
#include "TTree.h"

//#include "home/zmpur/HistoManager/HistoManager.h"
//#include "~/HistoManager.h"

Float_t DEGRAD=0.017453293;


struct PHYSDATA {float e, theta, phi;};
struct SABREDATA {float detectorIndex, theta, phi, energy, localx, localy; int ring, wedge;};

struct DetectorHit {
	int ring, wedge;
	double energy;
	//double ering, ewedge;
};

std::pair<double,double> getReconstructedAngles(int detectorIndex, int ring, int wedge, std::map<std::pair<int,int>,std::pair<double,double>> map);
std::map<std::pair<int,int>,std::pair<double,double>> readAngleMaps();
void readSingleAngleMap(ifstream& infile, std::map<std::pair<int,int>,std::pair<double,double>>& map);
void UpdateHistoAxes(HistoManager* histoman);

//std::map<std::pair<int,int>, std::pair<double,double>> sabre0_thetaphimap, sabre1_thetaphimap, sabre2_thetaphimap, sabre3_thetaphimap, sabre4_thetaphimap;

void analyze2BodyDetectorStepOutput(const char* input_filename, const char* output_rootfilename, const char* ntpname = "kin2"){
	
	std::map<std::pair<int,int>, std::pair<double,double>> sabre_thetaphimap = readAngleMaps();

	ifstream infile(input_filename);
	if(!infile.is_open()){
		cerr << "Error: Could not open file " << input_filename << endl;
		return;
	}

	TFile* outfile = new TFile(output_rootfilename,"RECREATE");
	HistoManager *histoman = new HistoManager(outfile);
	histoman->loadHistoConfig("RutherfordScattering_2body.HMConfig");

	TTree* kin2 = new TTree(ntpname, "3-body simulation tree");

	PHYSDATA physdata1, physdata2;//, physdata3, physdata4;		//holds the physics data for ejectile and recoil
	SABREDATA sabredata1, sabredata2;							//holds the sabre detection data for ejectile and recoil if detected by SABRE

	kin2->Branch("physdata1",&physdata1, "e/F:theta/F:phi/F");//branch to hold physics data for 
	kin2->Branch("physdata2",&physdata2, "e/F:theta/F:phi/F");

	kin2->Branch("sabredata1",&sabredata1,"detectorIndex/I:ring/I:wedge/I:theta/F:phi/F:energy/F");
	kin2->Branch("sabredata2",&sabredata2,"detectorIndex/I:ring/I:wedge/I:theta/F:phi/F:energy/F");

	//kin2->Branch();

	string line;
	//map<int,pair<int,int>> detector_hits;//ex: {404, {3,6}} means particle 4 in detector 4 in (r,w) = (3,6).
	map<int, DetectorHit> detector_hits;//ex: detector_hits[402] = {3,6,4.56} means that particle 4 (first digit) was detected in SABRE2 (last digit) with energy 4.56 in (ring=3,wedge=6)
	bool hit1 = false, hit2 = false;

	while(getline(infile,line)){

		//new:
		if(line == "-1"){//end event
			kin2->Fill();

			/*fill histograms w/ histomanager here*/
			//basic kinematics histograms:
			histoman->getHisto1D("hELab1")->Fill(physdata1.e);
			histoman->getHisto1D("hELab2")->Fill(physdata2.e);
			histoman->getHisto1D("hThetaLab1")->Fill(physdata1.theta);
			histoman->getHisto1D("hThetaLab2")->Fill(physdata2.theta);
			histoman->getHisto1D("hPhiLab1")->Fill(physdata1.phi);
			histoman->getHisto1D("hPhiLab2")->Fill(physdata2.phi);
			histoman->getHisto2D("hELabThetaLab_1")->Fill(physdata1.theta,physdata1.e);
			histoman->getHisto2D("hELabThetaLab_2")->Fill(physdata2.theta,physdata2.e);
			histoman->getHisto2D("hELabPhiLab_1")->Fill(physdata1.phi,physdata1.e);
			histoman->getHisto2D("hELabPhiLab_2")->Fill(physdata2.phi,physdata2.e);
			histoman->getHisto2D("hThetaLabPhiLab_1")->Fill(physdata1.theta,physdata1.phi);
			histoman->getHisto2D("hThetaLabPhiLab_2")->Fill(physdata2.theta,physdata2.phi);

			//SABRE ring/wedge hit summary histograms:
			if(sabredata1.detectorIndex == 0){
				histoman->getHisto1D("hSABRE0_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE0_WedgeHit")->Fill(sabredata1.wedge);
				histoman->getHisto1D("hSABRE0_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE0_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.energy);
				histoman->getHisto2D("hSABRE0_ESummaryRings")->Fill(sabredata1.ring, sabredata1.energy);
				histoman->getHisto2D("hSABRE0_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
			} else if(sabredata1.detectorIndex == 1) {
				histoman->getHisto1D("hSABRE1_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE1_WedgeHit")->Fill(sabredata1.wedge);
				histoman->getHisto1D("hSABRE1_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE1_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.energy);
				histoman->getHisto2D("hSABRE1_ESummaryRings")->Fill(sabredata1.ring, sabredata1.energy);
				histoman->getHisto2D("hSABRE1_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
			} else if(sabredata1.detectorIndex == 2){
				histoman->getHisto1D("hSABRE2_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE2_WedgeHit")->Fill(sabredata1.wedge);
				histoman->getHisto1D("hSABRE2_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE2_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.energy);
				histoman->getHisto2D("hSABRE2_ESummaryRings")->Fill(sabredata1.ring, sabredata1.energy);
				histoman->getHisto2D("hSABRE2_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
			} else if(sabredata1.detectorIndex == 3){
				histoman->getHisto1D("hSABRE3_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE3_WedgeHit")->Fill(sabredata1.wedge);
				histoman->getHisto1D("hSABRE3_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE3_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.energy);
				histoman->getHisto2D("hSABRE3_ESummaryRings")->Fill(sabredata1.ring, sabredata1.energy);
				histoman->getHisto2D("hSABRE3_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
			} else if(sabredata1.detectorIndex == 4){
				histoman->getHisto1D("hSABRE4_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE4_WedgeHit")->Fill(sabredata1.wedge);
				histoman->getHisto1D("hSABRE4_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE4_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.energy);
				histoman->getHisto2D("hSABRE4_ESummaryRings")->Fill(sabredata1.ring, sabredata1.energy);
				histoman->getHisto2D("hSABRE4_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
			}

			//SABRE ring/wedge hit summary histograms:
			if(sabredata2.detectorIndex == 0){
				histoman->getHisto1D("hSABRE0_RingHit")->Fill(sabredata2.ring);
				histoman->getHisto1D("hSABRE0_WedgeHit")->Fill(sabredata2.wedge);
				histoman->getHisto1D("hSABRE0_ESummary")->Fill(sabredata2.energy);
				histoman->getHisto2D("hSABRE0_ESummaryWedges")->Fill(sabredata2.wedge, sabredata2.energy);
				histoman->getHisto2D("hSABRE0_ESummaryRings")->Fill(sabredata2.ring, sabredata2.energy);
				histoman->getHisto2D("hSABRE0_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata2.theta, sabredata2.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
			} else if(sabredata2.detectorIndex == 1) {
				histoman->getHisto1D("hSABRE1_RingHit")->Fill(sabredata2.ring);
				histoman->getHisto1D("hSABRE1_WedgeHit")->Fill(sabredata2.wedge);
				histoman->getHisto1D("hSABRE1_ESummary")->Fill(sabredata2.energy);
				histoman->getHisto2D("hSABRE1_ESummaryWedges")->Fill(sabredata2.wedge, sabredata2.energy);
				histoman->getHisto2D("hSABRE1_ESummaryRings")->Fill(sabredata2.ring, sabredata2.energy);
				histoman->getHisto2D("hSABRE1_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata2.theta, sabredata2.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
			} else if(sabredata2.detectorIndex == 2){
				histoman->getHisto1D("hSABRE2_RingHit")->Fill(sabredata2.ring);
				histoman->getHisto1D("hSABRE2_WedgeHit")->Fill(sabredata2.wedge);
				histoman->getHisto1D("hSABRE2_ESummary")->Fill(sabredata2.energy);
				histoman->getHisto2D("hSABRE2_ESummaryWedges")->Fill(sabredata2.wedge, sabredata2.energy);
				histoman->getHisto2D("hSABRE2_ESummaryRings")->Fill(sabredata2.ring, sabredata2.energy);
				histoman->getHisto2D("hSABRE2_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata2.theta, sabredata2.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
			} else if(sabredata2.detectorIndex == 3){
				histoman->getHisto1D("hSABRE3_RingHit")->Fill(sabredata2.ring);
				histoman->getHisto1D("hSABRE3_WedgeHit")->Fill(sabredata2.wedge);
				histoman->getHisto1D("hSABRE3_ESummary")->Fill(sabredata2.energy);
				histoman->getHisto2D("hSABRE3_ESummaryWedges")->Fill(sabredata2.wedge, sabredata2.energy);
				histoman->getHisto2D("hSABRE3_ESummaryRings")->Fill(sabredata2.ring, sabredata2.energy);
				histoman->getHisto2D("hSABRE3_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata2.theta, sabredata2.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
			} else if(sabredata2.detectorIndex == 4){
				histoman->getHisto1D("hSABRE4_RingHit")->Fill(sabredata2.ring);
				histoman->getHisto1D("hSABRE4_WedgeHit")->Fill(sabredata2.wedge);
				histoman->getHisto1D("hSABRE4_ESummary")->Fill(sabredata2.energy);
				histoman->getHisto2D("hSABRE4_ESummaryWedges")->Fill(sabredata2.wedge, sabredata2.energy);
				histoman->getHisto2D("hSABRE4_ESummaryRings")->Fill(sabredata2.ring, sabredata2.energy);
				histoman->getHisto2D("hSABRE4_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata2.theta, sabredata2.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata2.localx, -sabredata2.localy);
			}



			detector_hits.clear();
			hit1 = false;
			hit2 = false;
			sabredata1 = {-1,-999,-999,-999,-999,-999};
			sabredata2 = {-1,-999,-999,-999,-999,-999};
		} else if(line.substr(0,2) == "10"){//particle1/ejectile detection
			stringstream ss(line);
			int detector_particle_id, ring, wedge;
			double energy, x, y;
			if(ss >> detector_particle_id >> ring >> wedge >> energy >> x >> y){
				int detector_index = detector_particle_id % 10;
				DetectorHit h;
				h.ring = ring;
				h.wedge = wedge;
				h.energy = energy;
				detector_hits[detector_particle_id] = h;
				//double energy = physdata1.e;
				std::pair<double,double> angles = getReconstructedAngles(detector_index, ring, wedge, sabre_thetaphimap);
				double theta_recon = angles.first;
				double phi_recon = angles.second;

				SABREDATA sabre_hit;
				sabre_hit.detectorIndex = detector_index;
				sabre_hit.ring = ring;
				sabre_hit.wedge = wedge;
				sabre_hit.theta = theta_recon;
				sabre_hit.phi = phi_recon;
				sabre_hit.energy = energy;
				sabre_hit.localx = x;
				sabre_hit.localy = y;

				if(!hit1 && !hit2){
					sabredata1 = sabre_hit;
					hit1 = true;
				}else if(hit1 && !hit2){
					sabredata2 = sabre_hit;
					hit2 = true;
				}

			}
		} else if(line.substr(0,2) == "20"){//particle2/recoil detection
			stringstream ss(line);
			int detector_particle_id, ring, wedge;
			double energy, x, y;
			if(ss >> detector_particle_id >> ring >> wedge >> energy >> x >> y){
				int detector_index = detector_particle_id % 10;
				DetectorHit h;
				h.ring = ring;
				h.wedge = wedge;
				h.energy = energy;
				detector_hits[detector_particle_id] = h;
				std::pair<double,double> angles = getReconstructedAngles(detector_index, ring, wedge, sabre_thetaphimap);
				double theta_recon = angles.first;
				double phi_recon = angles.second;

				SABREDATA sabre_hit;
				sabre_hit.detectorIndex = detector_index;
				sabre_hit.ring = ring;
				sabre_hit.wedge = wedge;
				sabre_hit.theta = theta_recon;
				sabre_hit.phi = phi_recon;
				sabre_hit.energy = energy;
				sabre_hit.localx = x;
				sabre_hit.localy = y;

				if(!hit1 && !hit2){
					sabredata1 = sabre_hit;
					hit1 = true;
				} else if(hit1 && !hit2){
					sabredata2 = sabre_hit;
					hit2 = true;
				}
			}
		} else {//kinematics
			stringstream ss(line);
			double dummythetacm;
			ss >> physdata1.e >> physdata1.theta >> physdata1.phi >> dummythetacm >> physdata2.e >> physdata2.theta >> physdata2.phi;
			if(physdata1.theta < 0) physdata1.theta += 180.;
			if(physdata1.phi < 0) physdata1.phi += 360.;
			if(physdata2.theta < 0) physdata2.theta += 180.;
			if(physdata2.phi < 0) physdata2.phi += 360.;
		}

	}

	outfile->cd();
	kin2->Write();

	/*Update Histogram Axes Here*/
	UpdateHistoAxes(histoman);

	histoman->WriteAll(true);
}

std::pair<double,double> getReconstructedAngles(int detectorIndex, int ring, int wedge, std::map<std::pair<int,int>,std::pair<double,double>> map){
	static const std::pair<int, int> offsets[] = {
		{112,40},	//detector0
		{96,32},	//detector1
		{80,16},	//detector2
		{64,24},	//detector3
		{48,0}		//detector4
	};

	if(detectorIndex < 0 || detectorIndex >= 5){
		return {-666.666,-666.666};
	}

	const auto& [ringOffset,wedgeOffset] = offsets[detectorIndex];
	return map[{ring+ringOffset,wedge+wedgeOffset}];
}

std::map<std::pair<int,int>,std::pair<double,double>> readAngleMaps(){
	const vector<string> filenames = {
		"SABRE0_phi306_anglemap.txt",
		"SABRE1_phi18_anglemap.txt",
		"SABRE2_phi234_anglemap.txt",
		"SABRE3_phi162_anglemap.txt",
		"SABRE4_phi90_anglemap.txt"
	};

	std::map<std::pair<int,int>,std::pair<double,double>> retmap;

	for(const auto& filename : filenames){
		ifstream infile(filename);
		if(!infile.is_open()){
			cerr << "Error: Failed to open file " << filename << endl;
			continue;
		}
		readSingleAngleMap(infile, retmap);
		infile.close();
	}

	return retmap;
}

void readSingleAngleMap(ifstream& infile, std::map<std::pair<int,int>,std::pair<double,double>>& map){
	string header;
	getline(infile,header);

	int ringChannel, wedgeChannel;
	double theta, phi;
	while(infile >> ringChannel >> wedgeChannel >> theta >> phi){
		map[{ringChannel,wedgeChannel}] = {theta,phi};
	}
}

void UpdateHistoAxes(HistoManager* histoman){
	histoman->getHisto2D("hSABRE0_hitsMapLocal")->GetXaxis()->SetTitle("+x <--------------------------------------------------> -x");
	histoman->getHisto2D("hSABRE0_hitsMapLocal")->GetYaxis()->SetTitle("+y <--------------------------------------------------> -y");
	histoman->getHisto2D("hSABRE0_hitsMapLocal")->GetXaxis()->CenterTitle();
	histoman->getHisto2D("hSABRE0_hitsMapLocal")->GetYaxis()->CenterTitle();

	histoman->getHisto2D("hSABRE1_hitsMapLocal")->GetXaxis()->SetTitle("+x <--------------------------------------------------> -x");
	histoman->getHisto2D("hSABRE1_hitsMapLocal")->GetYaxis()->SetTitle("+y <--------------------------------------------------> -y");
	histoman->getHisto2D("hSABRE1_hitsMapLocal")->GetXaxis()->CenterTitle();
	histoman->getHisto2D("hSABRE1_hitsMapLocal")->GetYaxis()->CenterTitle();

	histoman->getHisto2D("hSABRE2_hitsMapLocal")->GetXaxis()->SetTitle("+x <--------------------------------------------------> -x");
	histoman->getHisto2D("hSABRE2_hitsMapLocal")->GetYaxis()->SetTitle("+y <--------------------------------------------------> -y");
	histoman->getHisto2D("hSABRE2_hitsMapLocal")->GetXaxis()->CenterTitle();
	histoman->getHisto2D("hSABRE2_hitsMapLocal")->GetYaxis()->CenterTitle();

	histoman->getHisto2D("hSABRE3_hitsMapLocal")->GetXaxis()->SetTitle("+x <--------------------------------------------------> -x");
	histoman->getHisto2D("hSABRE3_hitsMapLocal")->GetYaxis()->SetTitle("+y <--------------------------------------------------> -y");
	histoman->getHisto2D("hSABRE3_hitsMapLocal")->GetXaxis()->CenterTitle();
	histoman->getHisto2D("hSABRE3_hitsMapLocal")->GetYaxis()->CenterTitle();

	histoman->getHisto2D("hSABRE4_hitsMapLocal")->GetXaxis()->SetTitle("+x <--------------------------------------------------> -x");
	histoman->getHisto2D("hSABRE4_hitsMapLocal")->GetYaxis()->SetTitle("+y <--------------------------------------------------> -y");
	histoman->getHisto2D("hSABRE4_hitsMapLocal")->GetXaxis()->CenterTitle();
	histoman->getHisto2D("hSABRE4_hitsMapLocal")->GetYaxis()->CenterTitle();

	histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->GetXaxis()->SetTitle("+x <--------------------------------------------------> -x");
	histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->GetYaxis()->SetTitle("+y <--------------------------------------------------> -y");
	histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->GetXaxis()->CenterTitle();
	histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->GetYaxis()->CenterTitle();

}
/*

*/
using namespace std;
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <vector>

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
#include "TVector3.h"
#include "TRotation.h"
#include "TRandom3.h"


Float_t DEGRAD=0.017453293;

int numRings = 16;
int numWedges = 8;


struct PHYSDATA {double e, theta, phi;};
struct SABREDATA {int detectorIndex=-666; double theta, phi, ringEnergy, wedgeEnergy, localx, localy; int ring, wedge;};


struct Corner{
	double x, y;
};

struct BinData{
	int ring, wege;
	std::vector<Corner> corners;
};

struct dp4Nucleus{
	//data to be passed to IMMMA_Tool_4 to set reaction nuclei
	int A;
	TString sym;
	double massMeV;

	void SetA(int newA){
		A = newA;
	}

	void SetSym(TString newsym){
		sym = newsym;
	}

	void SetMassMeV(double newmassMeV){
		massMeV = newmassMeV;
	}

	void SetAll(int newA, TString newsym, double newmassMeV){
		SetA(newA);
		SetSym(newsym);
		SetMassMeV(newmassMeV);
	}

};

struct Reaction4{
	dp4Nucleus beam;
	dp4Nucleus target;
	dp4Nucleus ejectile;
	dp4Nucleus recoil;
	dp4Nucleus breakup1;
	dp4Nucleus daughter1;
	dp4Nucleus breakup2;
	dp4Nucleus breakup3;

	double beamEnergy = 0.;
	double recoilExE = 0.;

	TString ToString(){
		TString retval = Form("%d%s(%d%s,%d%s)%d%s, %d%s -> %d%s + %d%s, %d%s -> %d%s + %d%s\tat E=%f to ExE=%f", target.A, target.sym.Data(), 
																												  beam.A, beam.sym.Data(),
																												  ejectile.A, ejectile.sym.Data(),
																												  recoil.A, recoil.sym.Data(),
																												  recoil.A, recoil.sym.Data(),
																												  daughter1.A, daughter1.sym.Data(),
																												  breakup1.A, breakup1.sym.Data(),
																												  daughter1.A, daughter1.sym.Data(), 
																												  breakup2.A, breakup2.sym.Data(),
																												  breakup3.A, breakup3.sym.Data(), 
																												  beamEnergy,
																												  recoilExE);

		return retval;
	};

};

struct CMCMinMax{
	//holds min max for CMConstants
	double min = std::numeric_limits<double>::max();
	double max = std::numeric_limits<double>::lowest();
	bool initialized = false;

	TString ToString(){
		return Form("min = %f, max = %f",min,max);
	}
};

static const std::pair<int, int> offsets[] = {
	{112,40},	//detector0 {ringOffset,wedgeOffset} // rings span channels 112 - 127, wedges span channels 40 - 47
	{96,32},	//detector1 {ringOffset,wedgeOffset} // rings span channels 96 - 111, wedges span channels 32 - 39
	{80,16},	//detector2 {ringOffset,wedgeOffset} // rings span channels 80 - 95, wedges span channels 16 - 23
	{64,24},	//detector3 {ringOffset,wedgeOffset} // rings span channels 64 - 69, wedges span channels 24 - 31
	{48,0}		//detector4 {ringOffset,wedgeOffset} // rings span channels 48 - 63, wedges span channels 0 - 7
};

double PERC = 0.10;
double LOWER = 1-PERC, UPPER = 1+PERC;

std::pair<double,double> getReconstructedAngles(int detectorIndex, int ring, int wedge, std::map<std::pair<int,int>,std::pair<double,double>> map);
std::map<std::pair<int,int>,std::pair<double,double>> readAngleMaps();
void readSingleAngleMap(ifstream& infile, std::map<std::pair<int,int>,std::pair<double,double>>& map);
void UpdateHistoAxes(HistoManager* histoman);
void DeterminePolygons(HistoManager *histoman);

bool parsePhysData(const std::string& line, PHYSDATA& pd1, PHYSDATA& pd2, PHYSDATA& pd3, PHYSDATA& pd4){
	std::istringstream iss(line);
	(iss >> pd1.e >> pd1.theta >> pd1.phi >> pd2.e >> pd2.theta >> pd2.phi >> pd3.e >> pd3.theta >> pd3.phi >> pd4.e >> pd4.theta >> pd4.phi);
	return true;
}

bool parseSABREData(const std::string& line, SABREDATA& sd){
	std::istringstream iss(line);
	int index, ring, wedge;
	double ringe, wedgee, x, y;
	iss >> index >> ring >> wedge >> ringe >> wedgee >> x >> y;
	// (iss >> sd.detectorIndex >> sd.ring >> sd.wedge >> sd.ringEnergy >> sd.wedgeEnergy >> sd.localx >> sd.localy);
	sd.detectorIndex = index%100;
	sd.ringEnergy = ringe;
	sd.wedgeEnergy = wedgee;
	sd.localx = x;
	sd.localy = y;
	sd.ring = ring;
	sd.wedge = wedge;
	sd.theta = -666.;//update this after the fact with angle map
	sd.phi = -666.;//update this after the fact with angle map
	return true;
}

void fillKinHistos(HistoManager* histoman, PHYSDATA& physdata1, PHYSDATA& physdata2, PHYSDATA& physdata3, PHYSDATA& physdata4){
	//basic kinematics histograms:
	histoman->getHisto1D("hELab1")->Fill(physdata1.e);
	histoman->getHisto1D("hELab2")->Fill(physdata2.e);
	histoman->getHisto1D("hELab3")->Fill(physdata3.e);
	histoman->getHisto1D("hELab4")->Fill(physdata4.e);
	histoman->getHisto1D("hThetaLab1")->Fill(physdata1.theta);
	histoman->getHisto1D("hThetaLab2")->Fill(physdata2.theta);
	histoman->getHisto1D("hThetaLab3")->Fill(physdata3.theta);
	histoman->getHisto1D("hThetaLab4")->Fill(physdata4.theta);
	histoman->getHisto1D("hPhiLab1")->Fill(physdata1.phi);
	histoman->getHisto1D("hPhiLab2")->Fill(physdata2.phi);
	histoman->getHisto1D("hPhiLab3")->Fill(physdata3.phi);
	histoman->getHisto1D("hPhiLab4")->Fill(physdata4.phi);
	histoman->getHisto2D("hELabThetaLab_1")->Fill(physdata1.theta,physdata1.e);
	histoman->getHisto2D("hELabThetaLab_2")->Fill(physdata2.theta,physdata2.e);
	histoman->getHisto2D("hELabThetaLab_3")->Fill(physdata3.theta,physdata3.e);
	histoman->getHisto2D("hELabThetaLab_4")->Fill(physdata4.theta,physdata4.e);
	histoman->getHisto2D("hELabPhiLab_1")->Fill(physdata1.phi,physdata1.e);
	histoman->getHisto2D("hELabPhiLab_2")->Fill(physdata2.phi,physdata2.e);
	histoman->getHisto2D("hELabPhiLab_3")->Fill(physdata3.phi,physdata3.e);
	histoman->getHisto2D("hELabPhiLab_4")->Fill(physdata4.phi,physdata4.e);
	histoman->getHisto2D("hThetaLabPhiLab_1")->Fill(physdata1.theta,physdata1.phi);
	histoman->getHisto2D("hThetaLabPhiLab_2")->Fill(physdata2.theta,physdata2.phi);
	histoman->getHisto2D("hThetaLabPhiLab_3")->Fill(physdata3.theta,physdata3.phi);
	histoman->getHisto2D("hThetaLabPhiLab_4")->Fill(physdata4.theta,physdata4.phi);
	//cout << "fillKinHistos test" << endl;
	TVector3 lab1(sin(DEGRAD*physdata3.theta)*cos(DEGRAD*physdata3.phi), sin(DEGRAD*physdata3.theta)*sin(DEGRAD*physdata3.phi), cos(DEGRAD*physdata3.theta)), lab2(sin(DEGRAD*physdata4.theta)*cos(DEGRAD*physdata4.phi), sin(DEGRAD*physdata4.theta)*sin(DEGRAD*physdata4.phi), cos(DEGRAD*physdata4.theta));
	histoman->getHisto1D("hBreakupLabAngle")->Fill(lab1.Angle(lab2)/DEGRAD);
}


void fillSABREHistos(HistoManager* histoman, SABREDATA& sabredata1, PHYSDATA &physdata1){
			int ringoffset = offsets[sabredata1.detectorIndex].first;
			int wedgeoffset = offsets[sabredata1.detectorIndex].second;
			int globalring = sabredata1.ring + ringoffset;
			int globalwedge = sabredata1.wedge + wedgeoffset;

			//rings vs wedges histogram:
			int zeroToFortyWedge = sabredata1.detectorIndex*numWedges + sabredata1.wedge;//0-7 for SABRE0, 8-15 for SABRE1, etc.
			histoman->getHisto2D("hSABRE_RingsVSWedges")->Fill(zeroToFortyWedge,sabredata1.ring);

			//hSABRE0_pixel_r112w40_ESummary
			//pixel histo:
			TString pixelhistoname = Form("hSABRE%d_pixel_r%dw%d_ESummary", sabredata1.detectorIndex, globalring, globalwedge);
			histoman->getHisto1D(pixelhistoname)->Fill(sabredata1.ringEnergy);

			//SABRE ring/wedge hit summary histograms:
			if(sabredata1.detectorIndex == 0){
				histoman->getHisto1D("hSABRE0_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE0_WedgeHit")->Fill(sabredata1.wedge);
				//histoman->getHisto1D("hSABRE0_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE0_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE0_ESummaryRings")->Fill(sabredata1.ring, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE0_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				//cout << "theta = " << sabredata1.theta << ", phi = " << sabredata1.phi << endl;
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				// histoman->getHisto2D("hSABRE0_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				// histoman->getHisto2D("hSABRE0_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				// histoman->getHisto2D("hSABRE0_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				// histoman->getHisto2D("hSABRE0_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
				// histoman->getHisto2D("hSABRE0_ringEDif")->Fill(sabredata1.ring, ringEDif);
				// histoman->getHisto2D("hSABRE0_wedgeEDif")->Fill(sabredata1.wedge, wedgeEDif);
				histoman->getHisto1D("hSABRE0_EDif")->Fill(sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto2DPoly("hSABRE_PixelMap")->Fill(sabredata1.localx, sabredata1.localy);
				histoman->getHisto3D("hSABRE0_PixelEDif")->Fill(sabredata1.wedge, sabredata1.ring, sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto3D("hSABRE0_ESummaryPixels")->Fill(sabredata1.wedge,sabredata1.ring,sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE0_ChannelHits")->Fill(globalring);
				histoman->getHisto1D("hSABRE0_ChannelHits")->Fill(globalwedge);
				histoman->getHisto1D("hSABRE0_ERingSummary")->Fill(sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE0_EWedgeSummary")->Fill(sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE0_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE0_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);
			} else if(sabredata1.detectorIndex == 1) {
				histoman->getHisto1D("hSABRE1_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE1_WedgeHit")->Fill(sabredata1.wedge);
				//histoman->getHisto1D("hSABRE1_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE1_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE1_ESummaryRings")->Fill(sabredata1.ring, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE1_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				// histoman->getHisto2D("hSABRE1_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				// histoman->getHisto2D("hSABRE1_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				// histoman->getHisto2D("hSABRE1_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				// histoman->getHisto2D("hSABRE1_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
				// histoman->getHisto2D("hSABRE1_ringEDif")->Fill(sabredata1.ring, ringEDif);
				// histoman->getHisto2D("hSABRE1_wedgeEDif")->Fill(sabredata1.wedge, wedgeEDif);
				histoman->getHisto1D("hSABRE1_EDif")->Fill(sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto2DPoly("hSABRE_PixelMap")->Fill(sabredata1.localx, sabredata1.localy);
				histoman->getHisto3D("hSABRE1_PixelEDif")->Fill(sabredata1.wedge, sabredata1.ring, sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto3D("hSABRE1_ESummaryPixels")->Fill(sabredata1.wedge,sabredata1.ring,sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE1_ChannelHits")->Fill(globalring);
				histoman->getHisto1D("hSABRE1_ChannelHits")->Fill(globalwedge);
				histoman->getHisto1D("hSABRE1_ERingSummary")->Fill(sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE1_EWedgeSummary")->Fill(sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE1_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE1_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);
			} else if(sabredata1.detectorIndex == 2){
				histoman->getHisto1D("hSABRE2_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE2_WedgeHit")->Fill(sabredata1.wedge);
				//histoman->getHisto1D("hSABRE2_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE2_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE2_ESummaryRings")->Fill(sabredata1.ring, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE2_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				// histoman->getHisto2D("hSABRE2_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				// histoman->getHisto2D("hSABRE2_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				// histoman->getHisto2D("hSABRE2_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				// histoman->getHisto2D("hSABRE2_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
				// histoman->getHisto2D("hSABRE2_ringEDif")->Fill(sabredata1.ring, ringEDif);
				// histoman->getHisto2D("hSABRE2_wedgeEDif")->Fill(sabredata1.wedge, wedgeEDif);
				histoman->getHisto1D("hSABRE2_EDif")->Fill(sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto2DPoly("hSABRE_PixelMap")->Fill(sabredata1.localx, sabredata1.localy);
				histoman->getHisto3D("hSABRE2_PixelEDif")->Fill(sabredata1.wedge, sabredata1.ring, sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto3D("hSABRE2_ESummaryPixels")->Fill(sabredata1.wedge,sabredata1.ring,sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE2_ChannelHits")->Fill(globalring);
				histoman->getHisto1D("hSABRE2_ChannelHits")->Fill(globalwedge);
				histoman->getHisto1D("hSABRE2_ERingSummary")->Fill(sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE2_EWedgeSummary")->Fill(sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE2_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE2_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);
			} else if(sabredata1.detectorIndex == 3){
				histoman->getHisto1D("hSABRE3_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE3_WedgeHit")->Fill(sabredata1.wedge);
				//histoman->getHisto1D("hSABRE3_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE3_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE3_ESummaryRings")->Fill(sabredata1.ring, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE3_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				// histoman->getHisto2D("hSABRE3_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				// histoman->getHisto2D("hSABRE3_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				// histoman->getHisto2D("hSABRE3_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				// histoman->getHisto2D("hSABRE3_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
				// histoman->getHisto2D("hSABRE3_ringEDif")->Fill(sabredata1.ring, ringEDif);
				// histoman->getHisto2D("hSABRE3_wedgeEDif")->Fill(sabredata1.wedge, wedgeEDif);
				histoman->getHisto1D("hSABRE3_EDif")->Fill(sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto2DPoly("hSABRE_PixelMap")->Fill(sabredata1.localx, sabredata1.localy);
				histoman->getHisto3D("hSABRE3_PixelEDif")->Fill(sabredata1.wedge, sabredata1.ring, sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto3D("hSABRE3_ESummaryPixels")->Fill(sabredata1.wedge,sabredata1.ring,sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE3_ChannelHits")->Fill(globalring);
				histoman->getHisto1D("hSABRE3_ChannelHits")->Fill(globalwedge);
				histoman->getHisto1D("hSABRE3_ERingSummary")->Fill(sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE3_EWedgeSummary")->Fill(sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE3_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE3_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);
			} else if(sabredata1.detectorIndex == 4){
				//update globalring value to account for reversed counting:
				//std::cout << "before: globalring = " << globalring << "\tlocal ring = " << sabredata1.ring << std::endl;
				//globalring = 63 - sabredata1.ring;
				//std::cout << "after: globalring = " << globalring << "\t local ring = " << sabredata1.ring << std::endl << std::endl;

				histoman->getHisto1D("hSABRE4_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE4_WedgeHit")->Fill(sabredata1.wedge);
				//histoman->getHisto1D("hSABRE4_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE4_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE4_ESummaryRings")->Fill(sabredata1.ring, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE4_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				// histoman->getHisto2D("hSABRE4_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				// histoman->getHisto2D("hSABRE4_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				// histoman->getHisto2D("hSABRE4_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				// histoman->getHisto2D("hSABRE4_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
				// histoman->getHisto2D("hSABRE4_ringEDif")->Fill(sabredata1.ring, ringEDif);
				// histoman->getHisto2D("hSABRE4_wedgeEDif")->Fill(sabredata1.wedge, wedgeEDif);
				histoman->getHisto1D("hSABRE4_EDif")->Fill(sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto2DPoly("hSABRE_PixelMap")->Fill(sabredata1.localx, sabredata1.localy);
				histoman->getHisto3D("hSABRE4_PixelEDif")->Fill(sabredata1.wedge, sabredata1.ring, sabredata1.ringEnergy-sabredata1.wedgeEnergy);
				histoman->getHisto3D("hSABRE4_ESummaryPixels")->Fill(sabredata1.wedge,sabredata1.ring,sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE4_ChannelHits")->Fill(globalring);
				histoman->getHisto1D("hSABRE4_ChannelHits")->Fill(globalwedge);
				histoman->getHisto1D("hSABRE4_ERingSummary")->Fill(sabredata1.ringEnergy);
				histoman->getHisto1D("hSABRE4_EWedgeSummary")->Fill(sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE4_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE4_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);
			}

			histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalring);
			histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalwedge);
			histoman->getHisto1D("hSABRE_RingChannelHits")->Fill(globalring);
			histoman->getHisto1D("hSABRE_WedgeChannelHits")->Fill(globalwedge);
			histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
			histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);

			//cout << "fillSABREHistos test" << endl;

			//all sabre histograms:
			// histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalring);
			// histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalwedge);
			// histoman->getHisto1D("hSABRE_RingChannelHits")->Fill(globalring);
			// histoman->getHisto1D("hSABRE_WedgeChannelHits")->Fill(globalwedge);
			// histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
			// histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);
}

std::vector<IMMMA_Tool_4*> prepareIMMMA_Tool_4s(std::vector<Reaction4>& reactions, bool output=false){
	if(reactions.size() == 0) std::cout << "No data in reactions" << std::endl;

	if(output) std::cout << std::endl;

	std::vector<IMMMA_Tool_4*> retvect;

	for(size_t i=0; i<reactions.size(); i++){
		
		IMMMA_Tool_4* tool = new IMMMA_Tool_4();

		tool->SetBeamNucleus(reactions[i].beam.A, reactions[i].beam.sym, reactions[i].beam.massMeV);
		tool->SetTargetNucleus(reactions[i].target.A, reactions[i].target.sym, reactions[i].target.massMeV);
		tool->SetEjectileNucleus(reactions[i].ejectile.A, reactions[i].ejectile.sym, reactions[i].ejectile.massMeV);
		tool->SetRecoilNucleus(reactions[i].recoil.A, reactions[i].recoil.sym, reactions[i].recoil.massMeV);
		tool->SetBreakup1Nucleus(reactions[i].breakup1.A, reactions[i].breakup1.sym, reactions[i].breakup1.massMeV);
		tool->SetDaughter1Nucleus(reactions[i].daughter1.A, reactions[i].daughter1.sym, reactions[i].daughter1.massMeV);
		tool->SetBreakup2Nucleus(reactions[i].breakup2.A, reactions[i].breakup2.sym, reactions[i].breakup2.massMeV);
		tool->SetBreakup3Nucleus(reactions[i].breakup3.A, reactions[i].breakup3.sym, reactions[i].breakup3.massMeV);
		

		tool->SetBeamEnergy(reactions[i].beamEnergy);
		tool->SetRecoilExE(reactions[i].recoilExE);

		tool->CalculateCMConstants();

		tool->SetVcm_bu1_boundsp(LOWER,UPPER);
		tool->SetVcm_bu2_boundsp(LOWER,UPPER);
		tool->SetVcm_bu3_boundsp(LOWER,UPPER);
		tool->SetKEcm_bu1_boundsp(LOWER,UPPER);
		tool->SetKEcm_bu2_boundsp(LOWER,UPPER);
		tool->SetKEcm_bu3_boundsp(LOWER,UPPER);
		tool->SetEcm1_boundsp(LOWER,UPPER);
		tool->SetEcm2_boundsp(LOWER,UPPER);

		if(output) std::cout << "Finished setting " << reactions[i].ToString() << std::endl;
		retvect.push_back(tool);
	}

	if(output) std::cout << std::endl;
	return retvect;
}

double getSPSEnergy(double kinEnMeV, double sigmaMeV=0.015){
	double smearedE;

	TRandom3* rand_gen = new TRandom3(0);
	smearedE = rand_gen->Gaus(kinEnMeV,sigmaMeV);

	delete rand_gen;
	return smearedE;
}

void UpdateCMCMinMax(CMCMinMax& minmax, double value){
	//checks if values < minmax.min or values > minmax.max and updates accordingly
	if(!minmax.initialized){
		minmax.min = value;
		minmax.max = value;
		minmax.initialized = true;
	} else {
		if(value < minmax.min){
			minmax.min = value;
		} else if(value > minmax.max){
			minmax.max = value;
		}
	}
}

double calculateSPS_ExE(double spsE, double spsTheta, double spsPhi, TMassTable& table){
	TLorentzVector beam, target, ejectile, recoil;
	double smearedSPSE = getSPSEnergy(spsE);
	double beamEnergy = 7.5;
	beam.SetPxPyPzE(0.,0.,sqrt(2*table.GetMassMeV("He",3)*beamEnergy),beamEnergy+table.GetMassMeV("He",3));
	target.SetPxPyPzE(0.,0.,0.,table.GetMassMeV("B",10));
	double pej = sqrt(2*smearedSPSE*table.GetMassMeV("He",4));
	ejectile.SetPxPyPzE(pej*sin(DEGRAD*spsTheta)*cos(DEGRAD*spsPhi), pej*sin(DEGRAD*spsTheta)*sin(DEGRAD*spsPhi), pej*cos(DEGRAD*spsTheta),smearedSPSE+table.GetMassMeV("He",4));
	recoil = beam + target - ejectile;

	return (recoil.M() - table.GetMassMeV("B",9));
}

void B10ha_3halfminus(const char* input_filename, const char* output_rootfilename, const char* ntpname = "kin4"){
	//this function assumes ejectile theta and phi are restricted to the SPS

	std::map<std::pair<int,int>,std::pair<double,double>> sabre_thetaphimap = readAngleMaps();

	ifstream infile(input_filename);
	if(!infile.is_open()){
		std::cerr << "Error: could not open file " << input_filename << std::endl;
		return;
	}

	TFile *outfile = new TFile(output_rootfilename,"RECREATE");
	HistoManager *histoman = new HistoManager(outfile);
	histoman->loadHistoConfig("HMConfig/_4body.HMconfig");

	//mass table
	TMassTable fMassTable;
	fMassTable.Init("../../config/masstable.dat");

	//masses
	Double_t beamMass = fMassTable.GetMassMeV("He",3);
	Double_t targetMass = fMassTable.GetMassMeV("B",10);
	Double_t ejectileMass = fMassTable.GetMassMeV("He",4);
	Double_t recoilMass = fMassTable.GetMassMeV("B",9);
	Double_t bu1Mass = fMassTable.GetMassMeV("H",1);
	Double_t daughter1Mass = fMassTable.GetMassMeV("Be",8);
	Double_t bu2Mass = fMassTable.GetMassMeV("He",4);
	Double_t bu3Mass = fMassTable.GetMassMeV("He",4);

	//masses for 9B -> 4He + 5Li, 5Li -> p + 4He:
	// Double_t beamMass = fMassTable.GetMassMeV("He",3);
	// Double_t targetMass = fMassTable.GetMassMeV("B",10);
	// Double_t ejectileMass = fMassTable.GetMassMeV("He",4);
	// Double_t recoilMass = fMassTable.GetMassMeV("B",9);
	// Double_t bu1Mass = fMassTable.GetMassMeV("He",4);
	// Double_t daughter1Mass = fMassTable.GetMassMeV("Li",5);
	// Double_t bu2Mass = fMassTable.GetMassMeV("H",1);
	// Double_t bu3Mass = fMassTable.GetMassMeV("He",4);

	//reactions
	//reaction 10B(3He,4He)9B, 9B -> p + 8Be, 8Be -> 4He + 4He
	std::vector<Reaction4> reactions;
	Reaction4 r,r2;
	r.beam.SetAll(3,"He",beamMass);
	r.target.SetAll(10, "B", targetMass);
	r.ejectile.SetAll(4, "He", ejectileMass);
	r.recoil.SetAll(9, "B", recoilMass);
	r.breakup1.SetAll(1, "H", bu1Mass);
	r.daughter1.SetAll(8, "Be", daughter1Mass);
	r.breakup2.SetAll(4, "He", bu2Mass);
	r.breakup3.SetAll(4, "He", bu3Mass);
	r.beamEnergy = 7.5;
	r.recoilExE = 2.186;
	reactions.push_back(r);

	//reaction 10B(3He,4He)9B, 9B -> 4He + 5Li, 5Li -> p + 4He
	// r2.beam.SetAll(3,"He",beamMass);
	// r2.target.SetAll(10,"B",targetMass);
	// r2.ejectile.SetAll(4,"He",ejectileMass);
	// r2.recoil.SetAll(9,"B",recoilMass);
	// r2.breakup1.SetAll(4,"He",bu1Mass);
	// r2.daughter1.SetAll(5,"Li",daughter1Mass);
	// r2.breakup2.SetAll(1,"p",bu2Mass);
	// r2.breakup3.SetAll(4,"He",bu3Mass);
	// r2.beamEnergy = 7.5;
	// r2.recoilExE = 0.;
	// reactions.push_back(r2);

	std::vector<IMMMA_Tool_4*> tools = prepareIMMMA_Tool_4s(reactions,true);
	std::cout << "VCM1  = " << tools[0]->GetExpectedVcm_bu1() << std::endl;
	std::cout << "VCM2  = " << tools[0]->GetExpectedVcm_bu2() << std::endl;
	std::cout << "VCM3  = " << tools[0]->GetExpectedVcm_bu3() << std::endl;
	std::cout << "KECM1 = " << tools[0]->GetExpectedKEcm_bu1() << std::endl;
	std::cout << "KECM2 = " << tools[0]->GetExpectedKEcm_bu2() << std::endl;
	std::cout << "KECM3 = " << tools[0]->GetExpectedKEcm_bu3() << std::endl;
	std::cout << "ECM1   = " << tools[0]->GetExpectedEcm1() << std::endl;
	std::cout << "ECM2   = " << tools[0]->GetExpectedEcm2() << std::endl;
	std::cout << endl;

	// std::cout << "9B -> 4He + 5Li, 5Li -> p + 4He:" << std::endl;
	// std::cout << "VCM1   = " << tools[1]->GetExpectedVcm_bu1() << std::endl;
	// std::cout << "VCM2   = " << tools[1]->GetExpectedVcm_bu2() << std::endl;
	// std::cout << "VCM3   = " << tools[1]->GetExpectedVcm_bu3() << std::endl;
	// std::cout << "KECM1  = " << tools[1]->GetExpectedKEcm_bu1() << std::endl;
	// std::cout << "KECM2  = " << tools[1]->GetExpectedKEcm_bu2() << std::endl;
	// std::cout << "KECM3  = " << tools[1]->GetExpectedKEcm_bu3() << std::endl;
	// std::cout << "ECM1    = " << tools[1]->GetExpectedEcm1() << std::endl;
	// std::cout << "ECM2    = " << tools[1]->GetExpectedEcm2() << std::endl;
	// std::cout << std::endl;

	TTree* kin4 = new TTree(ntpname, "10B 3par 9B_3halfminus");

	string line;

	std::cout << "Processing " << input_filename << "..." << std::endl;
	int count = 0;

	std::vector<std::string> eventLines;
	while(std::getline(infile,line)){
		if(line == "-1"){
			//eoev marker!
			PHYSDATA pd1, pd2, pd3, pd4;
			SABREDATA sd1, sd2, sd3;

			if(eventLines.size() == 1){
				
				//kinematics line only
				parsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				fillKinHistos(histoman,pd1,pd2,pd3,pd4);

			} else if(eventLines.size() == 2){
				
				parsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				fillKinHistos(histoman,pd1,pd2,pd3,pd4);

				//1 sabre hit (N-2)
				parseSABREData(eventLines[1],sd1);
				std::pair<double, double> anglepair = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}];
				sd1.theta = anglepair.first;
				sd1.phi = anglepair.second;

				fillSABREHistos(histoman,sd1,pd1);

				Double_t exe = calculateSPS_ExE(pd1.e, pd1.theta, pd1.phi, fMassTable);
				histoman->getHisto2D("hSABRE_SabreRingESumVsB9ExE")->Fill(exe, sd1.ringEnergy);
				histoman->getHisto1D("hSABRE_SabreRingESum")->Fill(sd1.ringEnergy);
				histoman->getHisto1D("hSPS_ExE")->Fill(exe);

				histoman->getHisto2D("h1par_SabreRingEVsB9ExE")->Fill(exe, sd1.ringEnergy);
				

			} else if(eventLines.size() == 3){

				parsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				fillKinHistos(histoman,pd1,pd2,pd3,pd4);
				
				//2 sabre hit (N-1)
				parseSABREData(eventLines[1],sd1);
				std::pair<double, double> anglepair = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}];
				sd1.theta = anglepair.first;
				sd1.phi = anglepair.second;

				fillSABREHistos(histoman,sd1,pd1);

				parseSABREData(eventLines[2],sd2);
				anglepair = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first, sd2.wedge+offsets[sd2.detectorIndex].second}];
				sd2.theta = anglepair.first;
				sd2.phi = anglepair.second;

				fillSABREHistos(histoman,sd2,pd2);

				Double_t exe = calculateSPS_ExE(pd1.e, pd1.theta, pd1.phi, fMassTable);
				histoman->getHisto2D("hSABRE_SabreRingESumVsB9ExE")->Fill(exe, sd1.ringEnergy+sd2.ringEnergy);
				histoman->getHisto1D("hSABRE_SabreRingESum")->Fill(sd1.ringEnergy+sd2.ringEnergy);
				histoman->getHisto1D("hSPS_ExE")->Fill(exe);

				histoman->getHisto2D("h2and3par_SabreRingEVsB9ExE")->Fill(exe, sd1.ringEnergy+sd2.ringEnergy);
				histoman->getHisto2D("h2par_SabreRingEVsB9ExE")->Fill(exe, sd1.ringEnergy+sd2.ringEnergy);
				

				//eventually do MMM here!
				std::array<CaseResult4,6> results = tools[0]->AnalyzeMMMEvent(pd1.e, pd1.theta, pd1.phi, sd1.ringEnergy, sd1.theta, sd1.phi, sd2.ringEnergy, sd2.theta, sd2.phi);
				
				for(size_t i=0; i<results.size(); i++){
					histoman->getHisto1D("hMMM_9BExE_2par_Be8_He4He4")->Fill(results[i].recoilExE);
				}


			} else if(eventLines.size() == 4){

				parsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				fillKinHistos(histoman,pd1,pd2,pd3,pd4);
				
				//3 sabre hit (N)
				parseSABREData(eventLines[1],sd1);
				std::pair<double,double> anglepair = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}];
				sd1.theta = anglepair.first;
				sd1.phi = anglepair.second;

				fillSABREHistos(histoman,sd1,pd1);

				parseSABREData(eventLines[2],sd2);
				anglepair = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first, sd2.wedge+offsets[sd2.detectorIndex].second}];
				sd2.theta = anglepair.first;
				sd2.phi = anglepair.second;

				fillSABREHistos(histoman,sd2,pd2);

				parseSABREData(eventLines[3],sd3);
				anglepair = sabre_thetaphimap[{sd3.ring+offsets[sd3.detectorIndex].first, sd3.wedge+offsets[sd3.detectorIndex].second}];
				sd3.theta = anglepair.first;
				sd3.phi = anglepair.second;

				fillSABREHistos(histoman,sd3,pd3);

				Double_t exe = calculateSPS_ExE(pd1.e, pd1.theta, pd1.phi, fMassTable);
				histoman->getHisto2D("hSABRE_SabreRingESumVsB9ExE")->Fill(exe, sd1.ringEnergy+sd2.ringEnergy+sd3.ringEnergy);
				histoman->getHisto1D("hSABRE_SabreRingESum")->Fill(sd1.ringEnergy+sd2.ringEnergy+sd3.ringEnergy);
				histoman->getHisto1D("hSPS_ExE")->Fill(exe);

				histoman->getHisto2D("h2and3par_SabreRingEVsB9ExE")->Fill(exe, sd1.ringEnergy+sd2.ringEnergy+sd3.ringEnergy);
				histoman->getHisto2D("h3par_SabreRingEVsB9ExE")->Fill(exe, sd1.ringEnergy+sd2.ringEnergy+sd3.ringEnergy);

				//eventually do IMM here!
				std::array<CaseResult4,6> results = tools[0]->AnalyzeIMMEvent(pd1.e, pd1.theta, pd1.phi, sd1.ringEnergy, sd1.theta, sd1.phi, sd2.ringEnergy, sd2.theta, sd2.phi, sd3.ringEnergy, sd3.theta, sd3.phi);


				for(size_t i=0; i<results.size(); i++){
					histoman->getHisto1D("hIMM_9BExE_3par_Be8_He4He4")->Fill(results[i].recoilExE);
				}



			} else if(eventLines.size() == 5){
				//4 sabre hit -- this should not happen if you properly use SPS angles for kin4mc ejectiles!
			} else {
				std::cerr << "Warning! Unexpected eventLines.size() = " << eventLines.size() << std::endl;
			}

			eventLines.clear();
			count += 1;
			if(count%100000==0) std::cout << "Processed " << count << " events..." << std::endl;

		} else {
			eventLines.push_back(line);
		}
	}

	outfile->cd();
	kin4->Write();

	UpdateHistoAxes(histoman);

	histoman->WriteAll(true);
	std::cout << "\nProcessed " << count << " events.\nROOT file saved to " << output_rootfilename << "\n";
}

std::pair<double,double> getReconstructedAngles(int detectorIndex, int ring, int wedge, std::map<std::pair<int,int>,std::pair<double,double>> map){

	if(detectorIndex < 0 || detectorIndex >= 5){
		return {-666.666,-666.666};
	}

	const auto& [ringOffset,wedgeOffset] = offsets[detectorIndex];
	return map[{ring+ringOffset,wedge+wedgeOffset}];
}

std::map<std::pair<int,int>,std::pair<double,double>> readAngleMaps(){
	const vector<string> filenames = {
		"../../anglemaps/SABRE0_phi306_anglemap.txt",
		"../../anglemaps/SABRE1_phi18_anglemap.txt",
		"../../anglemaps/SABRE2_phi234_anglemap.txt",
		"../../anglemaps/SABRE3_phi162_anglemap.txt",
		"../../anglemaps/SABRE4_phi90_anglemap.txt"
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

std::pair<int,int> findClosestRingWedge(Double_t theta_in, Double_t phi_in, const std::map<std::pair<int, int>, std::pair<double, double>>& angleMap){
	double minDist = std::numeric_limits<double>::max();
	std::pair<int,int> bestRingWedge = {-1,-1};

	for(const auto& entry : angleMap){
		const auto& [ring,wedge] = entry.first;
		const auto& [theta,phi] = entry.second;

		double dTheta = theta - theta_in;
		double dPhi = phi - phi_in;
		double distSq = dTheta*dTheta + dPhi*dPhi;

		if(distSq < minDist){
			minDist = distSq;
			bestRingWedge = {ring,wedge};
		}
	}

	return bestRingWedge;

}

void readSingleAngleMap(ifstream& infile, std::map<std::pair<int,int>,std::pair<double,double>>& map){
	string header;
	getline(infile,header);

	int ringChannel, wedgeChannel;
	double theta, phi;
	while(infile >> ringChannel >> wedgeChannel >> theta >> phi){
		map[{ringChannel,wedgeChannel}] = {theta,phi};
		//cout << "theta = " << theta << ", phi = " << phi << endl;
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

	// vector<TString> histonames = {
	// 	"hSABRE0_ringEVSkinE",
	// 	"hSABRE1_ringEVSkinE",
	// 	"hSABRE2_ringEVSkinE",
	// 	"hSABRE3_ringEVSkinE",
	// 	"hSABRE4_ringEVSkinE",
	// 	"hSABRE0_wedgeEVSkinE",
	// 	"hSABRE1_wedgeEVSkinE",
	// 	"hSABRE2_wedgeEVSkinE",
	// 	"hSABRE3_wedgeEVSkinE",
	// 	"hSABRE4_wedgeEVSkinE",
	// 	"hSABRE0_ringThetaVSkinTheta",
	// 	"hSABRE1_ringThetaVSkinTheta",
	// 	"hSABRE2_ringThetaVSkinTheta",
	// 	"hSABRE3_ringThetaVSkinTheta",
	// 	"hSABRE4_ringThetaVSkinTheta",
	// 	"hSABRE0_wedgePhiVSkinPhi",
	// 	"hSABRE1_wedgePhiVSkinPhi",
	// 	"hSABRE2_wedgePhiVSkinPhi",
	// 	"hSABRE3_wedgePhiVSkinPhi",
	// 	"hSABRE4_wedgePhiVSkinPhi"
	// };

	// for(const auto& name : histonames){
	// 	histoman->getHisto2D(name)->GetXaxis()->SetTitle("SABRESim");
	// 	histoman->getHisto2D(name)->GetYaxis()->SetTitle("Kin2mc");
	// 	histoman->getHisto2D(name)->GetXaxis()->CenterTitle();
	// 	histoman->getHisto2D(name)->GetYaxis()->CenterTitle();

	// 	//histoman->getHisto2D(name)->SetOption("SQUARE");
	// }

	vector<TString> histonames = {
		"hSABRE0_EDif",
		"hSABRE1_EDif",
		"hSABRE2_EDif",
		"hSABRE3_EDif",
		"hSABRE4_EDif"
	};

	for(const auto& name : histonames){
		histoman->getHisto1D(name)->GetXaxis()->SetTitle("E_{ring} - E_{wedge}, any pixel");
		histoman->getHisto1D(name)->GetXaxis()->CenterTitle();
	}

	histonames = {
		"hSABRE0_PixelEDif",
		"hSABRE1_PixelEDif",
		"hSABRE2_PixelEDif",
		"hSABRE3_PixelEDif",
		"hSABRE4_PixelEDif"
	};

	for(const auto& name : histonames){
		histoman->getHisto3D(name)->GetXaxis()->SetTitle("Ring Index");
		histoman->getHisto3D(name)->GetYaxis()->SetTitle("Wedge Index");
		histoman->getHisto3D(name)->GetZaxis()->SetTitle("E_{ring} - E_{wedge} for pixel (ring,wedge)");
		histoman->getHisto3D(name)->GetXaxis()->CenterTitle();
		histoman->getHisto3D(name)->GetYaxis()->CenterTitle();
		histoman->getHisto3D(name)->GetZaxis()->CenterTitle();
	}

	histonames = {
		"hSABRE0_ESummaryPixels",
		"hSABRE1_ESummaryPixels",
		"hSABRE2_ESummaryPixels",
		"hSABRE3_ESummaryPixels",
		"hSABRE4_ESummaryPixels"
	};

	for(const auto& name : histonames){
		histoman->getHisto3D(name)->GetXaxis()->SetTitle("Wedge Index");
		histoman->getHisto3D(name)->GetYaxis()->SetTitle("Ring Index");
		histoman->getHisto3D(name)->GetZaxis()->SetTitle("ESummary By Pixel (using ERing)");
		histoman->getHisto3D(name)->GetXaxis()->CenterTitle();
		histoman->getHisto3D(name)->GetYaxis()->CenterTitle();
		histoman->getHisto3D(name)->GetZaxis()->CenterTitle();
	}

}


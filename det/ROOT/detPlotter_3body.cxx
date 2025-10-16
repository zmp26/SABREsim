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

int numrings = 16;
int numwedges = 8;


struct PHYSDATA {double e, theta, phi;};
struct SABREDATA {int detectorIndex=-666; double theta, phi, ringEnergy, wedgeEnergy, localx, localy; int ring, wedge;};


struct Corner{
	double x, y;
};

struct BinData{
	int ring, wege;
	std::vector<Corner> corners;
};

struct dp3Nucleus{
	//data to be passed to IMMMA_Tool_3 to set reaction nuclei
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

struct Reaction{
	dp3Nucleus beam;
	dp3Nucleus target;
	dp3Nucleus ejectile;
	dp3Nucleus recoil;
	dp3Nucleus breakup1;
	dp3Nucleus breakup2;

	double beamEnergy=0.;
	double recoilExE=0.;

	TString ToString(){
		TString retval = Form("%d%s(%d%s,%d%s)%d%s at E=%f to ExE=%f",target.A,target.sym.Data(),beam.A,beam.sym.Data(),ejectile.A,ejectile.sym.Data(),recoil.A,recoil.sym.Data(),beamEnergy,recoilExE);
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

				histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalring);
				histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalwedge);
				histoman->getHisto1D("hSABRE_RingChannelHits")->Fill(globalring);
				histoman->getHisto1D("hSABRE_WedgeChannelHits")->Fill(globalwedge);
				histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);
			}
			//cout << "fillSABREHistos test" << endl;

			//all sabre histograms:
			// histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalring);
			// histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalwedge);
			// histoman->getHisto1D("hSABRE_RingChannelHits")->Fill(globalring);
			// histoman->getHisto1D("hSABRE_WedgeChannelHits")->Fill(globalwedge);
			// histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
			// histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);
}

std::vector<IMMMA_Tool_3*> prepareIMMMA_Tool_3s(std::vector<Reaction>& reactions, bool output=false){
	if(reactions.size() == 0) cout << "No data in reactions" << endl;

	if(output) cout << endl;
	// TMassTable fMassTable;
	// fMassTable.Init("/mnt/e/kinematics/IMMMA_Tool/threebody/masstable.dat");

	std::vector<IMMMA_Tool_3*> retvect;

	for(size_t i=0; i<reactions.size(); i++){
		IMMMA_Tool_3* tool = new IMMMA_Tool_3();
		tool->SetBeamNucleus(reactions[i].beam.A,reactions[i].beam.sym, reactions[i].beam.massMeV);
		tool->SetTargetNucleus(reactions[i].target.A,reactions[i].target.sym, reactions[i].target.massMeV);
		tool->SetEjectileNucleus(reactions[i].ejectile.A,reactions[i].ejectile.sym, reactions[i].ejectile.massMeV);
		tool->SetRecoilNucleus(reactions[i].recoil.A,reactions[i].recoil.sym, reactions[i].recoil.massMeV);
		tool->SetBreakup1Nucleus(reactions[i].breakup1.A,reactions[i].breakup1.sym, reactions[i].breakup1.massMeV);
		tool->SetBreakup2Nucleus(reactions[i].breakup2.A,reactions[i].breakup2.sym, reactions[i].breakup2.massMeV);

		tool->SetBeamEnergy(reactions[i].beamEnergy);
		tool->SetRecoilExE(reactions[i].recoilExE);

		tool->CalculateCMConstants();
		//tool->SetMeanToExpected_All();

		tool->SetVcm_bu1_boundsp(LOWER,UPPER);
		tool->SetVcm_bu2_boundsp(LOWER,UPPER);
		tool->SetKEcm_bu1_boundsp(LOWER,UPPER);
		tool->SetKEcm_bu2_boundsp(LOWER,UPPER);
		tool->SetEcm_boundsp(LOWER,UPPER);

		if(output) cout << "Finished setting " << reactions[i].ToString() << endl;
		retvect.push_back(tool);
	}

	if(output) cout << endl;
	return retvect;
}

double getSPSEnergy(double kinEnMeV, double sigmaMeV=0.015){
	double smearedE;

	TRandom3 *rand_gen = new TRandom3(0);
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
	target.SetPxPyPzE(0.,0.,0.,table.GetMassMeV("Li",7));
	double pej = sqrt(2*smearedSPSE*table.GetMassMeV("He",4));
	ejectile.SetPxPyPzE(pej*sin(DEGRAD*spsTheta)*cos(DEGRAD*spsPhi), pej*sin(DEGRAD*spsTheta)*sin(DEGRAD*spsPhi), pej*cos(DEGRAD*spsTheta),smearedSPSE+table.GetMassMeV("He",4));
	recoil = beam + target - ejectile;

	return (recoil.M() - table.GetMassMeV("Li",6));
}

void LiFha_3plus(const char* input_filename, const char* output_rootfilename, const char* ntpname = "kin3"){
	//this function assumes ejectile theta and phi are restricted to the SPS
	std::map<std::pair<int,int>,std::pair<double,double>> sabre_thetaphimap = readAngleMaps();

	ifstream infile(input_filename);
	if(!infile.is_open()){
		cerr << "Error: Could not open file " << input_filename << endl;
		return;
	}

	TFile *outfile = new TFile(output_rootfilename,"RECREATE");
	HistoManager *histoman = new HistoManager(outfile);
	histoman->loadHistoConfig("HMConfig/LiFha_2par_3plus.HMConfig");

	//mass table
	TMassTable fMassTable;
	//fMassTable.Init("/mnt/e/kinematics/IMMMA_Tool/threebody/masstable.dat");//uncomment this for surface laptop!
	//fMassTable.Init("/home/zmpur/IMMMA_Tool/threebody/masstable.dat"); //uncomment this for DESKTOP!
	fMassTable.Init("../../config/masstable.dat");


	Double_t beamMass = fMassTable.GetMassMeV("He",3);
	Double_t targetMass = fMassTable.GetMassMeV("Li",7);
	Double_t ejectileMass = fMassTable.GetMassMeV("He",4);
	Double_t recoilMass = fMassTable.GetMassMeV("Li",6);
	Double_t bu1Mass = fMassTable.GetMassMeV("He",4);
	Double_t bu2Mass = fMassTable.GetMassMeV("H",2);
	//reactions:
	std::vector<Reaction> reactions;
	Reaction r;
	r.beam.SetAll(3,"He",beamMass);
	r.target.SetAll(7,"Li",targetMass);
	r.ejectile.SetAll(4,"He",ejectileMass);
	r.recoil.SetAll(6,"Li",recoilMass);
	r.breakup1.SetAll(4,"He",bu1Mass);
	r.breakup2.SetAll(2,"H",bu2Mass);
	r.beamEnergy = 7.5;
	r.recoilExE = 2.186;
	reactions.push_back(r);

	std::vector<IMMMA_Tool_3*> tools = prepareIMMMA_Tool_3s(reactions,true);
	cout << "VCM1 (4He)  = " << tools[0]->GetExpectedVcm_bu1() << endl;
	cout << "VCM2 (2H)   = " << tools[0]->GetExpectedVcm_bu2() << endl;
	cout << "KECM1 (4He) = " << tools[0]->GetExpectedKEcm_bu1() << endl;
	cout << "KECM2 (2H)  = " << tools[0]->GetExpectedKEcm_bu2() << endl;
	cout << "ECM         = " << tools[0]->GetExpectedEcm() << endl;
	cout << endl;

	TTree* kin3 = new TTree(ntpname,"LiFha 2par 6Li_3plus");

	// PHYSDATA physdata1, physdata2, physdata3, physdata4;
	// SABREDATA sabredata1, sabredata2;

	string line;
	bool hit1 = false, hit2 = false;

	cout << "Processing " << input_filename << "..." << endl;
	int count = 0;

	std::vector<std::string> eventLines;
	while(std::getline(infile,line)){
		if(line == "-1"){
			PHYSDATA pd1, pd2, pd3, pd4;
			SABREDATA sd1, sd2;
			if(eventLines.size()==1){
				//only kinematics line
				parsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				fillKinHistos(histoman,pd1,pd2,pd3,pd4);
				//kin IMM to get thetacm and phicm quickly:
				std::pair<CaseResult,CaseResult> kin = tools[0]->AnalyzeEventIMM(pd1.e,pd1.theta,pd1.phi,pd3.e,pd3.theta,pd3.phi,pd4.e,pd4.theta,pd4.phi);
				histoman->getHisto1D("hThetaCM_3")->Fill(kin.first.ThetaCM1);
				histoman->getHisto1D("hThetaCM_4")->Fill(kin.first.ThetaCM2);
				histoman->getHisto1D("hPhiCM_3")->Fill(kin.first.PhiCM1);
				histoman->getHisto1D("hPhiCM_4")->Fill(kin.first.PhiCM2);
				histoman->getHisto2D("hThetaCM3_vs_ThetaCM4")->Fill(kin.first.ThetaCM1,kin.first.ThetaCM2);
				histoman->getHisto2D("hCosThetaCM3_vs_CosThetaCM4")->Fill(TMath::Cos(kin.first.ThetaCM1*DEGRAD),TMath::Cos(kin.first.ThetaCM2*DEGRAD));
			} else if(eventLines.size()==2){
				//1 sabre hit
				parsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				parseSABREData(eventLines[1],sd1);
				sd1.theta = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].first;//wow this is ugly but it works
				sd1.phi = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].second;//wow this is ugly but it works
				fillKinHistos(histoman,pd1,pd2,pd3,pd4);
				fillSABREHistos(histoman,sd1,pd1);
				Double_t exe = calculateSPS_ExE(pd1.e,pd1.theta,pd1.phi,fMassTable);
				histoman->getHisto2D("hSABRE_SabreRingESumVsLi6ExE")->Fill(exe,sd1.ringEnergy);
				histoman->getHisto1D("hSABRE_SabreRingESum")->Fill(sd1.ringEnergy);
				histoman->getHisto1D("hSPS_ExE")->Fill(exe);
				//kin IMM to get thetacm and phicm quickly:
				std::pair<CaseResult,CaseResult> kin = tools[0]->AnalyzeEventIMM(pd1.e,pd1.theta,pd1.phi,pd3.e,pd3.theta,pd3.phi,pd4.e,pd4.theta,pd4.phi);
				histoman->getHisto1D("hThetaCM_3")->Fill(kin.first.ThetaCM1);
				histoman->getHisto1D("hThetaCM_4")->Fill(kin.first.ThetaCM2);
				histoman->getHisto1D("hPhiCM_3")->Fill(kin.first.PhiCM1);
				histoman->getHisto1D("hPhiCM_4")->Fill(kin.first.PhiCM2);
				histoman->getHisto2D("hThetaCM3_vs_ThetaCM4")->Fill(kin.first.ThetaCM1,kin.first.ThetaCM2);
				histoman->getHisto2D("hCosThetaCM3_vs_CosThetaCM4")->Fill(TMath::Cos(kin.first.ThetaCM1*DEGRAD),TMath::Cos(kin.first.ThetaCM2*DEGRAD));

				//do MMM here
				std::pair<CaseResult,CaseResult> results = tools[0]->AnalyzeEventMMM(pd1.e, pd1.theta, pd1.phi, sd1.ringEnergy, sd1.theta, sd1.phi);
				CaseResult case1 = results.first;
				CaseResult case2 = results.second;
				//fill angle between lab vector and VCM (in lab frame) vector histograms here
				histoman->getHisto1D("h1par_RecMissMassExE")->Fill(case1.recInvMass-recoilMass);
				histoman->getHisto1D("h1par_RecMissMassExE")->Fill(case2.recInvMass-recoilMass);
				histoman->getHisto1D("h1par_LabVCMAngle1")->Fill(case1.breakup1_LabAngleWRTVCM);
				histoman->getHisto1D("h1par_LabVCMAngle1")->Fill(case2.breakup1_LabAngleWRTVCM);
				histoman->getHisto1D("h1par_LabVCMAngle2")->Fill(case1.breakup2_LabAngleWRTVCM);
				histoman->getHisto1D("h1par_LabVCMAngle2")->Fill(case2.breakup2_LabAngleWRTVCM);
				histoman->getHisto1D("h1par_LabVCMAngleSum")->Fill(case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM);
				histoman->getHisto1D("h1par_LabVCMAngleSum")->Fill(case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM);
				histoman->getHisto1D("h1par_CosLabVCMAngleSum")->Fill(TMath::Cos((case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM)*DEGRAD));
				histoman->getHisto1D("h1par_CosLabVCMAngleSum")->Fill(TMath::Cos((case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM)*DEGRAD));
				histoman->getHisto1D("h1par_CosLabVCMAngle1")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h1par_CosLabVCMAngle1")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h1par_CosLabVCMAngle2")->Fill(TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h1par_CosLabVCMAngle2")->Fill(TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto2D("h1par_LabVCMAngle1vs2")->Fill(case1.breakup1_LabAngleWRTVCM,case1.breakup2_LabAngleWRTVCM);
				histoman->getHisto2D("h1par_LabVCMAngle1vs2")->Fill(case2.breakup1_LabAngleWRTVCM,case2.breakup2_LabAngleWRTVCM);
				histoman->getHisto2D("h1par_CosLabVsCosLab")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD),TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto2D("h1par_CosLabVsCosLab")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD),TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
				//fill angle between CM vector and VCM (in lab frame) vector histograms here
				histoman->getHisto1D("h1par_CMVCMAngle1")->Fill(case1.breakup1_CMAngleWRTVCM);
				histoman->getHisto1D("h1par_CMVCMAngle1")->Fill(case2.breakup1_CMAngleWRTVCM);
				histoman->getHisto1D("h1par_CMVCMAngle2")->Fill(case1.breakup2_CMAngleWRTVCM);
				histoman->getHisto1D("h1par_CMVCMAngle2")->Fill(case2.breakup2_CMAngleWRTVCM);
				histoman->getHisto1D("h1par_CMVCMAngleSum")->Fill(case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM);
				histoman->getHisto1D("h1par_CMVCMAngleSum")->Fill(case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM);
				histoman->getHisto1D("h1par_CosCMVCMAngleSum")->Fill(TMath::Cos((case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM)*DEGRAD));
				histoman->getHisto1D("h1par_CosCMVCMAngleSum")->Fill(TMath::Cos((case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM)*DEGRAD));
				histoman->getHisto1D("h1par_CosCMVCMAngle1")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h1par_CosCMVCMAngle1")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h1par_CosCMVCMAngle2")->Fill(TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h1par_CosCMVCMAngle2")->Fill(TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));
				histoman->getHisto2D("h1par_CMVCMAngle1vs2")->Fill(case1.breakup1_CMAngleWRTVCM,case1.breakup2_CMAngleWRTVCM);
				histoman->getHisto2D("h1par_CMVCMAngle1vs2")->Fill(case2.breakup1_CMAngleWRTVCM,case2.breakup2_CMAngleWRTVCM);
				histoman->getHisto2D("h1par_CosCMVsCosCM")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM),TMath::Cos(case1.breakup2_CMAngleWRTVCM));
				histoman->getHisto2D("h1par_CosCMVsCosCM")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM),TMath::Cos(case2.breakup2_CMAngleWRTVCM));

				histoman->getHisto2D("h1par_SabreRingEVsLi6ExE")->Fill(exe,sd1.ringEnergy);
				histoman->getHisto2D("h1and2par_SabreRingEVsLi6ExE")->Fill(exe,sd1.ringEnergy);

			} else if(eventLines.size()==3){
				//2 sabre hit
				parsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				parseSABREData(eventLines[1],sd1);
				sd1.theta = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first,sd1.wedge+offsets[sd1.detectorIndex].second}].first;
				sd1.phi = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first,sd1.wedge+offsets[sd1.detectorIndex].second}].second;
				parseSABREData(eventLines[2],sd2);
				sd2.theta = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first,sd2.wedge+offsets[sd2.detectorIndex].second}].first;
				sd2.phi = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first,sd2.wedge+offsets[sd2.detectorIndex].second}].second;
				fillKinHistos(histoman,pd1,pd2,pd3,pd4);
				fillSABREHistos(histoman,sd1,pd3);
				fillSABREHistos(histoman,sd2,pd4);
				Double_t exe = calculateSPS_ExE(pd1.e,pd1.theta,pd1.phi,fMassTable);
				histoman->getHisto2D("hSABRE_SabreRingESumVsLi6ExE")->Fill(exe,sd1.ringEnergy+sd2.ringEnergy);
				histoman->getHisto1D("hSABRE_SabreRingESum")->Fill(sd1.ringEnergy+sd2.ringEnergy);
				histoman->getHisto1D("hSPS_ExE")->Fill(exe);
				//kin IMM to get thetacm and phicm quickly:
				std::pair<CaseResult,CaseResult> kin = tools[0]->AnalyzeEventIMM(pd1.e,pd1.theta,pd1.phi,pd3.e,pd3.theta,pd3.phi,pd4.e,pd4.theta,pd4.phi);
				histoman->getHisto1D("hThetaCM_3")->Fill(kin.first.ThetaCM1);
				histoman->getHisto1D("hThetaCM_4")->Fill(kin.first.ThetaCM2);
				histoman->getHisto1D("hPhiCM_3")->Fill(kin.first.PhiCM1);
				histoman->getHisto1D("hPhiCM_4")->Fill(kin.first.PhiCM2);
				histoman->getHisto2D("hThetaCM3_vs_ThetaCM4")->Fill(kin.first.ThetaCM1,kin.first.ThetaCM2);
				histoman->getHisto2D("hCosThetaCM3_vs_CosThetaCM4")->Fill(TMath::Cos(kin.first.ThetaCM1*DEGRAD),TMath::Cos(kin.first.ThetaCM2*DEGRAD));
				//do IMM here
				std::pair<CaseResult,CaseResult> results = tools[0]->AnalyzeEventIMM(pd1.e,pd1.theta,pd1.phi,sd1.ringEnergy,sd1.theta,sd1.phi,sd2.ringEnergy,sd2.theta,sd2.phi);
				CaseResult case1 = results.first;
				CaseResult case2 = results.second;
				//fill histograms here:
				histoman->getHisto1D("h2par_RecInvMassExE")->Fill(case1.recInvMass-recoilMass);
				histoman->getHisto1D("h2par_RecInvMassExE")->Fill(case2.recInvMass-recoilMass);
				histoman->getHisto1D("h2par_LabVCMAngle1")->Fill(case1.breakup1_LabAngleWRTVCM);
				histoman->getHisto1D("h2par_LabVCMAngle1")->Fill(case2.breakup1_LabAngleWRTVCM);
				histoman->getHisto1D("h2par_LabVCMAngle2")->Fill(case1.breakup2_LabAngleWRTVCM);
				histoman->getHisto1D("h2par_LabVCMAngle2")->Fill(case2.breakup2_LabAngleWRTVCM);
				histoman->getHisto1D("h2par_LabVCMAngleSum")->Fill(case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM);
				histoman->getHisto1D("h2par_LabVCMAngleSum")->Fill(case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM);
				histoman->getHisto1D("h2par_CosLabVCMAngleSum")->Fill(TMath::Cos((case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM)*DEGRAD));
				histoman->getHisto1D("h2par_CosLabVCMAngleSum")->Fill(TMath::Cos((case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM)*DEGRAD));
				histoman->getHisto1D("h2par_CosLabVCMAngle1")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h2par_CosLabVCMAngle1")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h2par_CosLabVCMAngle2")->Fill(TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h2par_CosLabVCMAngle2")->Fill(TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto2D("h2par_LabVCMAngle1vs2")->Fill(case1.breakup1_LabAngleWRTVCM,case1.breakup2_LabAngleWRTVCM);
				histoman->getHisto2D("h2par_LabVCMAngle1vs2")->Fill(case2.breakup1_LabAngleWRTVCM,case2.breakup2_LabAngleWRTVCM);
				histoman->getHisto2D("h2par_CosLabVsCosLab")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD),TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
				histoman->getHisto2D("h2par_CosLabVsCosLab")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD),TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
				//fill angle between CM vector and VCM (in lab frame) vector histograms here
				histoman->getHisto1D("h2par_CMVCMAngle1")->Fill(case1.breakup1_CMAngleWRTVCM);
				histoman->getHisto1D("h2par_CMVCMAngle1")->Fill(case2.breakup1_CMAngleWRTVCM);
				histoman->getHisto1D("h2par_CMVCMAngle2")->Fill(case1.breakup2_CMAngleWRTVCM);
				histoman->getHisto1D("h2par_CMVCMAngle2")->Fill(case2.breakup2_CMAngleWRTVCM);
				histoman->getHisto1D("h2par_CMVCMAngleSum")->Fill(case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM);
				histoman->getHisto1D("h2par_CMVCMAngleSum")->Fill(case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM);
				histoman->getHisto1D("h2par_CosCMVCMAngleSum")->Fill(TMath::Cos((case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM)*DEGRAD));
				histoman->getHisto1D("h2par_CosCMVCMAngleSum")->Fill(TMath::Cos((case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM)*DEGRAD));
				histoman->getHisto1D("h2par_CosCMVCMAngle1")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h2par_CosCMVCMAngle1")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h2par_CosCMVCMAngle2")->Fill(TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
				histoman->getHisto1D("h2par_CosCMVCMAngle2")->Fill(TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));
				histoman->getHisto2D("h2par_CMVCMAngle1vs2")->Fill(case1.breakup1_CMAngleWRTVCM,case1.breakup2_CMAngleWRTVCM);
				histoman->getHisto2D("h2par_CMVCMAngle1vs2")->Fill(case2.breakup1_CMAngleWRTVCM,case2.breakup2_CMAngleWRTVCM);
				histoman->getHisto2D("h2par_CosCMVsCosCM")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM),TMath::Cos(case1.breakup2_CMAngleWRTVCM));
				histoman->getHisto2D("h2par_CosCMVsCosCM")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM),TMath::Cos(case2.breakup2_CMAngleWRTVCM));

				histoman->getHisto2D("h2par_SabreRingESumVsLi6ExE_3plus")->Fill(exe,sd1.ringEnergy+sd2.ringEnergy);
				histoman->getHisto2D("h1and2par_SabreRingEVsLi6ExE")->Fill(exe,sd1.ringEnergy+sd2.ringEnergy);

				int ringoffset = offsets[sd1.detectorIndex].first;
				int wedgeoffset = offsets[sd1.detectorIndex].second;
				int globalring = sd1.ring + ringoffset;
				int globalwedge = sd1.wedge + wedgeoffset;

				//if(sd1.detectorIndex==4) globalring = 63 - sd1.ring;

				histoman->getHisto1D("h2SABRE_ChannelHits_3plus")->Fill(globalring);
				histoman->getHisto1D("h2SABRE_ChannelHits_3plus")->Fill(globalwedge);
				histoman->getHisto2D("h2SABRE_ChannelESummary_3plus")->Fill(globalring,sd1.ringEnergy);
				histoman->getHisto2D("h2SABRE_ChannelESummary_3plus")->Fill(globalwedge,sd1.wedgeEnergy);

				ringoffset = offsets[sd2.detectorIndex].first;
				wedgeoffset = offsets[sd2.detectorIndex].second;
				globalring = sd2.ring + ringoffset;
				globalwedge = sd2.wedge + wedgeoffset;

				histoman->getHisto1D("h2SABRE_ChannelHits_3plus")->Fill(globalring);
				histoman->getHisto1D("h2SABRE_ChannelHits_3plus")->Fill(globalwedge);
				histoman->getHisto2D("h2SABRE_ChannelESummary_3plus")->Fill(globalring,sd2.ringEnergy);
				histoman->getHisto2D("h2SABRE_ChannelESummary_3plus")->Fill(globalwedge,sd2.wedgeEnergy);

				histoman->getHisto2D("h2par_LabThetaVsLabPhi1")->Fill(case1.PhiLab1,case1.ThetaLab1);
				histoman->getHisto2D("h2par_LabThetaVsLabPhi1")->Fill(case2.PhiLab1,case2.ThetaLab1);
				histoman->getHisto2D("h2par_LabThetaVsLabPhi2")->Fill(case1.PhiLab2,case1.ThetaLab2);
				histoman->getHisto2D("h2par_LabThetaVsLabPhi2")->Fill(case2.PhiLab2,case2.ThetaLab2);
				histoman->getHisto2D("h2par_LabTheta1VsLabTheta2")->Fill(case1.ThetaLab2,case1.ThetaLab1);
				histoman->getHisto2D("h2par_LabTheta1VsLabTheta2")->Fill(case2.ThetaLab2,case2.ThetaLab1);
				histoman->getHisto2D("h2par_LabPhi1VsLabPhi2")->Fill(case1.PhiLab2,case1.PhiLab1);
				histoman->getHisto2D("h2par_LabPhi1VsLabPhi2")->Fill(case2.PhiLab2,case2.PhiLab1);
				histoman->getHisto2D("h2par_LabTheta1VsLabPhi2")->Fill(case1.PhiLab2,case1.ThetaLab1);
				histoman->getHisto2D("h2par_LabTheta1VsLabPhi2")->Fill(case2.PhiLab2,case2.ThetaLab1);
				histoman->getHisto2D("h2par_LabTheta2VsLabPhi1")->Fill(case1.PhiLab1,case1.ThetaLab2);
				histoman->getHisto2D("h2par_LabTheta2VsLabPhi1")->Fill(case2.PhiLab1,case2.ThetaLab2);

				histoman->getHisto2D("h2par_CMThetaVsCMPhi1")->Fill(case1.PhiCM1,case1.ThetaCM1);
				histoman->getHisto2D("h2par_CMThetaVsCMPhi1")->Fill(case2.PhiCM1,case2.ThetaCM1);
				histoman->getHisto2D("h2par_CMThetaVsCMPhi2")->Fill(case1.PhiCM2,case1.ThetaCM2);
				histoman->getHisto2D("h2par_CMThetaVsCMPhi2")->Fill(case2.PhiCM2,case2.ThetaCM2);
				histoman->getHisto2D("h2par_CMTheta1VsCMTheta2")->Fill(case1.ThetaCM2,case1.ThetaCM1);
				histoman->getHisto2D("h2par_CMTheta1VsCMTheta2")->Fill(case2.ThetaCM2,case2.ThetaCM1);
				histoman->getHisto2D("h2par_CMPhi1VsCMPhi2")->Fill(case1.PhiCM2,case1.PhiCM1);
				histoman->getHisto2D("h2par_CMPhi1VsCMPhi2")->Fill(case2.PhiCM2,case2.PhiCM1);
				histoman->getHisto2D("h2par_CMTheta1VsCMPhi2")->Fill(case1.PhiCM2,case1.ThetaCM1);
				histoman->getHisto2D("h2par_CMTheta1VsCMPhi2")->Fill(case2.PhiCM2,case2.ThetaCM1);
				histoman->getHisto2D("h2par_CMTheta2VsCMPhi1")->Fill(case1.PhiCM1,case1.ThetaCM2);
				histoman->getHisto2D("h2par_CMTheta2VsCMPhi1")->Fill(case2.PhiCM1,case2.ThetaCM2);

				//ExE Cut checks here:
				double checkval = case1.recInvMass-recoilMass;
				if(checkval >= 2.05 && checkval <= 2.35){
					histoman->getHisto1D("h2par_LabVCMAngle1_ExECut")->Fill(case1.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle1_ExECut")->Fill(case2.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_ExECut")->Fill(case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_ExECut")->Fill(case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_ExECut")->Fill(case1.breakup1_LabAngleWRTVCM,case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_ExECut")->Fill(case2.breakup1_LabAngleWRTVCM,case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_ExECut")->Fill(case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_ExECut")->Fill(case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_ExECut")->Fill(TMath::Cos((case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_ExECut")->Fill(TMath::Cos((case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_ExECut")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_ExECut")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_ExECut")->Fill(TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_ExECut")->Fill(TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_ExECut")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_ExECut")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					//fill angle between CM vector and VCM (in lab frame) vector histograms here
					histoman->getHisto1D("h2par_CMVCMAngle1_ExECut")->Fill(case1.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle1_ExECut")->Fill(case2.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_ExECut")->Fill(case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_ExECut")->Fill(case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_ExECut")->Fill(case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_ExECut")->Fill(case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_ExECut")->Fill(case1.breakup1_CMAngleWRTVCM,case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_ExECut")->Fill(case2.breakup1_CMAngleWRTVCM,case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_ExECut")->Fill(TMath::Cos((case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_ExECut")->Fill(TMath::Cos((case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_ExECut")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_ExECut")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_ExECut")->Fill(TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_ExECut")->Fill(TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_ExECut")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_ExECut")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));

					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_ExECut")->Fill(case1.PhiLab1,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_ExECut")->Fill(case2.PhiLab1,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_ExECut")->Fill(case1.PhiLab2,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_ExECut")->Fill(case2.PhiLab2,case2.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_ExECut")->Fill(case1.ThetaLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_ExECut")->Fill(case2.ThetaLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_ExECut")->Fill(case1.PhiLab2,case1.PhiLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_ExECut")->Fill(case2.PhiLab2,case2.PhiLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_ExECut")->Fill(case1.PhiLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_ExECut")->Fill(case2.PhiLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_ExECut")->Fill(case1.PhiLab1,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_ExECut")->Fill(case2.PhiLab1,case2.ThetaLab2);

					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_ExECut")->Fill(case1.PhiCM1,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_ExECut")->Fill(case2.PhiCM1,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_ExECut")->Fill(case1.PhiCM2,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_ExECut")->Fill(case2.PhiCM2,case2.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_ExECut")->Fill(case1.ThetaCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_ExECut")->Fill(case2.ThetaCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_ExECut")->Fill(case1.PhiCM2,case1.PhiCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_ExECut")->Fill(case2.PhiCM2,case2.PhiCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_ExECut")->Fill(case1.PhiCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_ExECut")->Fill(case2.PhiCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_ExECut")->Fill(case1.PhiCM1,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_ExECut")->Fill(case2.PhiCM1,case2.ThetaCM2);
				} else if(checkval < 2.05){
					//ExECutLEFT
					histoman->getHisto1D("h2par_LabVCMAngle1_ExECutLEFT")->Fill(case1.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle1_ExECutLEFT")->Fill(case2.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_ExECutLEFT")->Fill(case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_ExECutLEFT")->Fill(case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_ExECutLEFT")->Fill(case1.breakup1_LabAngleWRTVCM,case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_ExECutLEFT")->Fill(case2.breakup1_LabAngleWRTVCM,case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_ExECutLEFT")->Fill(case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_ExECutLEFT")->Fill(case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_ExECutLEFT")->Fill(TMath::Cos((case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_ExECutLEFT")->Fill(TMath::Cos((case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_ExECutLEFT")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_ExECutLEFT")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_ExECutLEFT")->Fill(TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_ExECutLEFT")->Fill(TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_ExECutLEFT")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_ExECutLEFT")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					//fill angle between CM vector and VCM (in lab frame) vector histograms here
					histoman->getHisto1D("h2par_CMVCMAngle1_ExECutLEFT")->Fill(case1.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle1_ExECutLEFT")->Fill(case2.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_ExECutLEFT")->Fill(case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_ExECutLEFT")->Fill(case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_ExECutLEFT")->Fill(case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_ExECutLEFT")->Fill(case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_ExECutLEFT")->Fill(case1.breakup1_CMAngleWRTVCM,case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_ExECutLEFT")->Fill(case2.breakup1_CMAngleWRTVCM,case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_ExECutLEFT")->Fill(TMath::Cos((case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_ExECutLEFT")->Fill(TMath::Cos((case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_ExECutLEFT")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_ExECutLEFT")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_ExECutLEFT")->Fill(TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_ExECutLEFT")->Fill(TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_ExECutLEFT")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_ExECutLEFT")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));
					//AntiExECut
					histoman->getHisto1D("h2par_LabVCMAngle1_AntiExECut")->Fill(case1.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle1_AntiExECut")->Fill(case2.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_AntiExECut")->Fill(case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_AntiExECut")->Fill(case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_AntiExECut")->Fill(case1.breakup1_LabAngleWRTVCM,case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_AntiExECut")->Fill(case2.breakup1_LabAngleWRTVCM,case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_AntiExECut")->Fill(case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_AntiExECut")->Fill(case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_AntiExECut")->Fill(TMath::Cos((case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_AntiExECut")->Fill(TMath::Cos((case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_AntiExECut")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_AntiExECut")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_AntiExECut")->Fill(TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_AntiExECut")->Fill(TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_AntiExECut")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_AntiExECut")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					//fill angle between CM vector and VCM (in lab frame) vector histograms here
					histoman->getHisto1D("h2par_CMVCMAngle1_AntiExECut")->Fill(case1.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle1_AntiExECut")->Fill(case2.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_AntiExECut")->Fill(case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_AntiExECut")->Fill(case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_AntiExECut")->Fill(case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_AntiExECut")->Fill(case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_AntiExECut")->Fill(case1.breakup1_CMAngleWRTVCM,case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_AntiExECut")->Fill(case2.breakup1_CMAngleWRTVCM,case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_AntiExECut")->Fill(TMath::Cos((case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_AntiExECut")->Fill(TMath::Cos((case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_AntiExECut")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_AntiExECut")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_AntiExECut")->Fill(TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_AntiExECut")->Fill(TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_AntiExECut")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_AntiExECut")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));

					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_ExECutLEFT")->Fill(case1.PhiLab1,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_ExECutLEFT")->Fill(case2.PhiLab1,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_ExECutLEFT")->Fill(case1.PhiLab2,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_ExECutLEFT")->Fill(case2.PhiLab2,case2.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_ExECutLEFT")->Fill(case1.ThetaLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_ExECutLEFT")->Fill(case2.ThetaLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_ExECutLEFT")->Fill(case1.PhiLab2,case1.PhiLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_ExECutLEFT")->Fill(case2.PhiLab2,case2.PhiLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_ExECutLEFT")->Fill(case1.PhiLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_ExECutLEFT")->Fill(case2.PhiLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_ExECutLEFT")->Fill(case1.PhiLab1,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_ExECutLEFT")->Fill(case2.PhiLab1,case2.ThetaLab2);

					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_AntiExECut")->Fill(case1.PhiLab1,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_AntiExECut")->Fill(case2.PhiLab1,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_AntiExECut")->Fill(case1.PhiLab2,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_AntiExECut")->Fill(case2.PhiLab2,case2.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_AntiExECut")->Fill(case1.ThetaLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_AntiExECut")->Fill(case2.ThetaLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_AntiExECut")->Fill(case1.PhiLab2,case1.PhiLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_AntiExECut")->Fill(case2.PhiLab2,case2.PhiLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_AntiExECut")->Fill(case1.PhiLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_AntiExECut")->Fill(case2.PhiLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_AntiExECut")->Fill(case1.PhiLab1,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_AntiExECut")->Fill(case2.PhiLab1,case2.ThetaLab2);

					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_ExECutLEFT")->Fill(case1.PhiCM1,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_ExECutLEFT")->Fill(case2.PhiCM1,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_ExECutLEFT")->Fill(case1.PhiCM2,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_ExECutLEFT")->Fill(case2.PhiCM2,case2.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_ExECutLEFT")->Fill(case1.ThetaCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_ExECutLEFT")->Fill(case2.ThetaCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_ExECutLEFT")->Fill(case1.PhiCM2,case1.PhiCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_ExECutLEFT")->Fill(case2.PhiCM2,case2.PhiCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_ExECutLEFT")->Fill(case1.PhiCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_ExECutLEFT")->Fill(case2.PhiCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_ExECutLEFT")->Fill(case1.PhiCM1,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_ExECutLEFT")->Fill(case2.PhiCM1,case2.ThetaCM2);

					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_AntiExECut")->Fill(case1.PhiCM1,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_AntiExECut")->Fill(case2.PhiCM1,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_AntiExECut")->Fill(case1.PhiCM2,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_AntiExECut")->Fill(case2.PhiCM2,case2.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_AntiExECut")->Fill(case1.ThetaCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_AntiExECut")->Fill(case2.ThetaCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_AntiExECut")->Fill(case1.PhiCM2,case1.PhiCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_AntiExECut")->Fill(case2.PhiCM2,case2.PhiCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_AntiExECut")->Fill(case1.PhiCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_AntiExECut")->Fill(case2.PhiCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_AntiExECut")->Fill(case1.PhiCM1,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_AntiExECut")->Fill(case2.PhiCM1,case2.ThetaCM2);
				} else if(checkval > 2.35){
					//ExECutRIGHT
					histoman->getHisto1D("h2par_LabVCMAngle1_ExECutRIGHT")->Fill(case1.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle1_ExECutRIGHT")->Fill(case2.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_ExECutRIGHT")->Fill(case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_ExECutRIGHT")->Fill(case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_ExECutRIGHT")->Fill(case1.breakup1_LabAngleWRTVCM,case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_ExECutRIGHT")->Fill(case2.breakup1_LabAngleWRTVCM,case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_ExECutRIGHT")->Fill(case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_ExECutRIGHT")->Fill(case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_ExECutRIGHT")->Fill(TMath::Cos((case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_ExECutRIGHT")->Fill(TMath::Cos((case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_ExECutRIGHT")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_ExECutRIGHT")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_ExECutRIGHT")->Fill(TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_ExECutRIGHT")->Fill(TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_ExECutRIGHT")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_ExECutRIGHT")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					//fill angle between CM vector and VCM (in lab frame) vector histograms here
					histoman->getHisto1D("h2par_CMVCMAngle1_ExECutRIGHT")->Fill(case1.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle1_ExECutRIGHT")->Fill(case2.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_ExECutRIGHT")->Fill(case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_ExECutRIGHT")->Fill(case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_ExECutRIGHT")->Fill(case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_ExECutRIGHT")->Fill(case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_ExECutRIGHT")->Fill(case1.breakup1_CMAngleWRTVCM,case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_ExECutRIGHT")->Fill(case2.breakup1_CMAngleWRTVCM,case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_ExECutRIGHT")->Fill(TMath::Cos((case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_ExECutRIGHT")->Fill(TMath::Cos((case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_ExECutRIGHT")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_ExECutRIGHT")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_ExECutRIGHT")->Fill(TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_ExECutRIGHT")->Fill(TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_ExECutRIGHT")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_ExECutRIGHT")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));
					//AntiExECut
					histoman->getHisto1D("h2par_LabVCMAngle1_AntiExECut")->Fill(case1.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle1_AntiExECut")->Fill(case2.breakup1_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_AntiExECut")->Fill(case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngle2_AntiExECut")->Fill(case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_AntiExECut")->Fill(case1.breakup1_LabAngleWRTVCM,case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto2D("h2par_LabVCMAngle1vs2_AntiExECut")->Fill(case2.breakup1_LabAngleWRTVCM,case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_AntiExECut")->Fill(case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_LabVCMAngleSum_AntiExECut")->Fill(case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM);
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_AntiExECut")->Fill(TMath::Cos((case1.breakup1_LabAngleWRTVCM+case1.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngleSum_AntiExECut")->Fill(TMath::Cos((case2.breakup1_LabAngleWRTVCM+case2.breakup2_LabAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_AntiExECut")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle1_AntiExECut")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_AntiExECut")->Fill(TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosLabVCMAngle2_AntiExECut")->Fill(TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_AntiExECut")->Fill(TMath::Cos(case1.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_LabAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosLabVsCosLab_AntiExECut")->Fill(TMath::Cos(case2.breakup1_LabAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_LabAngleWRTVCM*DEGRAD));
					//fill angle between CM vector and VCM (in lab frame) vector histograms here
					histoman->getHisto1D("h2par_CMVCMAngle1_AntiExECut")->Fill(case1.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle1_AntiExECut")->Fill(case2.breakup1_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_AntiExECut")->Fill(case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngle2_AntiExECut")->Fill(case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_AntiExECut")->Fill(case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CMVCMAngleSum_AntiExECut")->Fill(case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_AntiExECut")->Fill(case1.breakup1_CMAngleWRTVCM,case1.breakup2_CMAngleWRTVCM);
					histoman->getHisto2D("h2par_CMVCMAngle1vs2_AntiExECut")->Fill(case2.breakup1_CMAngleWRTVCM,case2.breakup2_CMAngleWRTVCM);
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_AntiExECut")->Fill(TMath::Cos((case1.breakup1_CMAngleWRTVCM+case1.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngleSum_AntiExECut")->Fill(TMath::Cos((case2.breakup1_CMAngleWRTVCM+case2.breakup2_CMAngleWRTVCM)*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_AntiExECut")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle1_AntiExECut")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_AntiExECut")->Fill(TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto1D("h2par_CosCMVCMAngle2_AntiExECut")->Fill(TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_AntiExECut")->Fill(TMath::Cos(case1.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case1.breakup2_CMAngleWRTVCM*DEGRAD));
					histoman->getHisto2D("h2par_CosCMVsCosCM_AntiExECut")->Fill(TMath::Cos(case2.breakup1_CMAngleWRTVCM*DEGRAD), TMath::Cos(case2.breakup2_CMAngleWRTVCM*DEGRAD));

					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_ExECutRIGHT")->Fill(case1.PhiLab1,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_ExECutRIGHT")->Fill(case2.PhiLab1,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_ExECutRIGHT")->Fill(case1.PhiLab2,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_ExECutRIGHT")->Fill(case2.PhiLab2,case2.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_ExECutRIGHT")->Fill(case1.ThetaLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_ExECutRIGHT")->Fill(case2.ThetaLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_ExECutRIGHT")->Fill(case1.PhiLab2,case1.PhiLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_ExECutRIGHT")->Fill(case2.PhiLab2,case2.PhiLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_ExECutRIGHT")->Fill(case1.PhiLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_ExECutRIGHT")->Fill(case2.PhiLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_ExECutRIGHT")->Fill(case1.PhiLab1,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_ExECutRIGHT")->Fill(case2.PhiLab1,case2.ThetaLab2);

					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_AntiExECut")->Fill(case1.PhiLab1,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi1_AntiExECut")->Fill(case2.PhiLab1,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_AntiExECut")->Fill(case1.PhiLab2,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabThetaVsLabPhi2_AntiExECut")->Fill(case2.PhiLab2,case2.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_AntiExECut")->Fill(case1.ThetaLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabTheta2_AntiExECut")->Fill(case2.ThetaLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_AntiExECut")->Fill(case1.PhiLab2,case1.PhiLab1);
					histoman->getHisto2D("h2par_LabPhi1VsLabPhi2_AntiExECut")->Fill(case2.PhiLab2,case2.PhiLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_AntiExECut")->Fill(case1.PhiLab2,case1.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta1VsLabPhi2_AntiExECut")->Fill(case2.PhiLab2,case2.ThetaLab1);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_AntiExECut")->Fill(case1.PhiLab1,case1.ThetaLab2);
					histoman->getHisto2D("h2par_LabTheta2VsLabPhi1_AntiExECut")->Fill(case2.PhiLab1,case2.ThetaLab2);

					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_ExECutRIGHT")->Fill(case1.PhiCM1,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_ExECutRIGHT")->Fill(case2.PhiCM1,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_ExECutRIGHT")->Fill(case1.PhiCM2,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_ExECutRIGHT")->Fill(case2.PhiCM2,case2.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_ExECutRIGHT")->Fill(case1.ThetaCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_ExECutRIGHT")->Fill(case2.ThetaCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_ExECutRIGHT")->Fill(case1.PhiCM2,case1.PhiCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_ExECutRIGHT")->Fill(case2.PhiCM2,case2.PhiCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_ExECutRIGHT")->Fill(case1.PhiCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_ExECutRIGHT")->Fill(case2.PhiCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_ExECutRIGHT")->Fill(case1.PhiCM1,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_ExECutRIGHT")->Fill(case2.PhiCM1,case2.ThetaCM2);

					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_AntiExECut")->Fill(case1.PhiCM1,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi1_AntiExECut")->Fill(case2.PhiCM1,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_AntiExECut")->Fill(case1.PhiCM2,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMThetaVsCMPhi2_AntiExECut")->Fill(case2.PhiCM2,case2.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_AntiExECut")->Fill(case1.ThetaCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMTheta2_AntiExECut")->Fill(case2.ThetaCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_AntiExECut")->Fill(case1.PhiCM2,case1.PhiCM1);
					histoman->getHisto2D("h2par_CMPhi1VsCMPhi2_AntiExECut")->Fill(case2.PhiCM2,case2.PhiCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_AntiExECut")->Fill(case1.PhiCM2,case1.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta1VsCMPhi2_AntiExECut")->Fill(case2.PhiCM2,case2.ThetaCM1);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_AntiExECut")->Fill(case1.PhiCM1,case1.ThetaCM2);
					histoman->getHisto2D("h2par_CMTheta2VsCMPhi1_AntiExECut")->Fill(case2.PhiCM1,case2.ThetaCM2);
				}

			} else if(eventLines.size()==4){
				//3 sabre hit, this should not happen with kin3mc if you restrict theta and phi to SPS!
			} else {
				cerr << "Warning: Unexpected eventLines.size() = " << eventLines.size() << endl;
			}

			eventLines.clear();
			count += 1;
			if(count%100000==0) cout << "Processed " << count << " events..." << endl;

		} else {
			eventLines.push_back(line);
		}
	}

	outfile->cd();
	kin3->Write();

	UpdateHistoAxes(histoman);

	histoman->WriteAll(true);
	cout << endl;
	cout << "Processed " << count << " events.\n";
	cout << "ROOT file saved to " << output_rootfilename << "\n";
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

void DeterminePolygons(HistoManager* histoman){
	gROOT->SetBatch(kTRUE);

	std::ifstream corners[5];
	for(int i=0; i<5; i++){
		TString fn = Form("corners/SABRE%d_corners.txt",i);
		corners[i].open(fn.Data());
	}

	TH2Poly* hpoly = histoman->getHisto2DPoly("hSABRE_PixelMap");

	for(int det=0; det<5; det++){
		int ring, wedge;
		double x[4], y[4];
		while(corners[det] >> ring >> wedge >> x[0] >> y[0] >> x[1] >> y[1] >> x[2] >> y[2] >> x[3] >> y[3]){
			hpoly->AddBin(4,x,y);
		}

	}

	hpoly->GetXaxis()->SetLimits(-0.2,0.2);
	hpoly->GetYaxis()->SetLimits(-0.2,0.2);


	for(int i=0; i<5; i++){
		corners[i].close();
	}

}
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
#include "TVector3.h"
#include "TRotation.h"

//#include "home/zmpur/HistoManager/HistoManager.h"
//#include "~/HistoManager.h"

Float_t DEGRAD=0.017453293;


struct PHYSDATA {double e, theta, phi;};
struct SABREDATA {int detectorIndex=-666; double theta, phi, ringEnergy, wedgeEnergy, localx, localy; int ring, wedge;};

// struct DetectorHit {
// 	int ring, wedge;
// 	double ringEnergy, wedgeEnergy;
// 	//double ering, ewedge;
// };

struct Corner{
	double x, y;
};

struct BinData{
	int ring, wege;
	std::vector<Corner> corners;
};

static const std::pair<int, int> offsets[] = {
	{112,40},	//detector0 {ringOffset,wedgeOffset}
	{96,32},	//detector1 {ringOffset,wedgeOffset}
	{80,16},	//detector2 {ringOffset,wedgeOffset}
	{64,24},	//detector3 {ringOffset,wedgeOffset}
	{48,0}		//detector4 {ringOffset,wedgeOffset}
};

double PERC = 0.1;
double LOWER = 1-PERC, UPPER = 1+PERC;
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

std::pair<double,double> getReconstructedAngles(int detectorIndex, int ring, int wedge, std::map<std::pair<int,int>,std::pair<double,double>> map);
std::map<std::pair<int,int>,std::pair<double,double>> readAngleMaps();
void readSingleAngleMap(ifstream& infile, std::map<std::pair<int,int>,std::pair<double,double>>& map);
void UpdateHistoAxes(HistoManager* histoman);
void DeterminePolygons(HistoManager *histoman);
bool parsePhysData(const std::string& line, PHYSDATA& pd1, PHYSDATA& pd2){
	std::istringstream iss(line);
	double thetacm;
	(iss >> pd1.e >> pd1.theta >> pd1.phi >> thetacm >> pd2.e >> pd2.theta >> pd2.phi);
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

void fillKinHistos(HistoManager* histoman, PHYSDATA& physdata1, PHYSDATA& physdata2){
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
	//cout << "fillKinHistos test" << endl;
}

void fillSABREHistos(HistoManager* histoman, SABREDATA& sabredata1, PHYSDATA &physdata1){
			int ringoffset = offsets[sabredata1.detectorIndex].first;
			int wedgeoffset = offsets[sabredata1.detectorIndex].second;
			int globalring = sabredata1.ring + ringoffset;
			int globalwedge = sabredata1.wedge + wedgeoffset;

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
				histoman->getHisto2D("hSABRE0_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE0_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE0_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				histoman->getHisto2D("hSABRE0_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
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
				if(sabredata1.ring==0) { histoman->getHisto1D("hSABRE0_ring0_E")->Fill(sabredata1.ringEnergy); histoman->getHisto1D("hSABRE_ring0_E")->Fill(sabredata1.ringEnergy);}
			} else if(sabredata1.detectorIndex == 1) {
				histoman->getHisto1D("hSABRE1_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE1_WedgeHit")->Fill(sabredata1.wedge);
				//histoman->getHisto1D("hSABRE1_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE1_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE1_ESummaryRings")->Fill(sabredata1.ring, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE1_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE1_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE1_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE1_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				histoman->getHisto2D("hSABRE1_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
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
				if(sabredata1.ring==0) { histoman->getHisto1D("hSABRE1_ring0_E")->Fill(sabredata1.ringEnergy); histoman->getHisto1D("hSABRE_ring0_E")->Fill(sabredata1.ringEnergy);}
			} else if(sabredata1.detectorIndex == 2){
				histoman->getHisto1D("hSABRE2_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE2_WedgeHit")->Fill(sabredata1.wedge);
				//histoman->getHisto1D("hSABRE2_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE2_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE2_ESummaryRings")->Fill(sabredata1.ring, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE2_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE2_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE2_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE2_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				histoman->getHisto2D("hSABRE2_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
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
				if(sabredata1.ring==0) { histoman->getHisto1D("hSABRE2_ring0_E")->Fill(sabredata1.ringEnergy); histoman->getHisto1D("hSABRE_ring0_E")->Fill(sabredata1.ringEnergy);}
			} else if(sabredata1.detectorIndex == 3){
				histoman->getHisto1D("hSABRE3_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE3_WedgeHit")->Fill(sabredata1.wedge);
				//histoman->getHisto1D("hSABRE3_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE3_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE3_ESummaryRings")->Fill(sabredata1.ring, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE3_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE3_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE3_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE3_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				histoman->getHisto2D("hSABRE3_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
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
				if(sabredata1.ring==0) { histoman->getHisto1D("hSABRE3_ring0_E")->Fill(sabredata1.ringEnergy); histoman->getHisto1D("hSABRE_ring0_E")->Fill(sabredata1.ringEnergy);}
				if(sabredata1.ring==7){ histoman->getHisto1D("hSABRE3_Ring7Summary")->Fill(sabredata1.ringEnergy);}
				if(sabredata1.ring==8){ histoman->getHisto1D("hSABRE3_Ring8Summary")->Fill(sabredata1.ringEnergy);}
				if(sabredata1.ring==9){ histoman->getHisto1D("hSABRE3_Ring9Summary")->Fill(sabredata1.ringEnergy);}
			} else if(sabredata1.detectorIndex == 4){
				histoman->getHisto1D("hSABRE4_RingHit")->Fill(sabredata1.ring);
				histoman->getHisto1D("hSABRE4_WedgeHit")->Fill(sabredata1.wedge);
				//histoman->getHisto1D("hSABRE4_ESummary")->Fill(sabredata1.energy);
				histoman->getHisto2D("hSABRE4_ESummaryWedges")->Fill(sabredata1.wedge, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE4_ESummaryRings")->Fill(sabredata1.ring, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE4_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE_AngleHitsMap")->Fill(sabredata1.theta, sabredata1.phi);
				histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sabredata1.localx, -sabredata1.localy);
				histoman->getHisto2D("hSABRE4_ringEVSkinE")->Fill(physdata1.e, sabredata1.ringEnergy);
				histoman->getHisto2D("hSABRE4_wedgeEVSkinE")->Fill(physdata1.e, sabredata1.wedgeEnergy);
				histoman->getHisto2D("hSABRE4_ringThetaVSkinTheta")->Fill(physdata1.theta, sabredata1.theta);
				histoman->getHisto2D("hSABRE4_wedgePhiVSkinPhi")->Fill(physdata1.phi, sabredata1.phi);
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
				if(sabredata1.ring==0) { histoman->getHisto1D("hSABRE4_ring0_E")->Fill(sabredata1.ringEnergy); histoman->getHisto1D("hSABRE_ring0_E")->Fill(sabredata1.ringEnergy);}
			}
			//cout << "fillSABREHistos test" << endl;

			//all sabre histograms:
			histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalring);
			histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalwedge);
			histoman->getHisto1D("hSABRE_RingChannelHits")->Fill(globalring);
			histoman->getHisto1D("hSABRE_WedgeChannelHits")->Fill(globalwedge);
			histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalring,sabredata1.ringEnergy);
			histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalwedge,sabredata1.wedgeEnergy);
} 


//std::map<std::pair<int,int>, std::pair<double,double>> sabre0_thetaphimap, sabre1_thetaphimap, sabre2_thetaphimap, sabre3_thetaphimap, sabre4_thetaphimap;

double getSPSEnergy(double kinEnMeV, double sigmaMeV=0.015){
	double smearedE;

	TRandom3 *rand_gen = new TRandom3(0);
	smearedE = rand_gen->Gaus(kinEnMeV,sigmaMeV);

	delete rand_gen;

	return smearedE;
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

void analyze2BodyDetectorStepOutput(const char* input_filename, const char* output_rootfilename, const char* ntpname = "kin2"){
	
	std::map<std::pair<int,int>, std::pair<double,double>> sabre_thetaphimap = readAngleMaps();

	ifstream infile(input_filename);
	if(!infile.is_open()){
		cerr << "Error: Could not open file " << input_filename << endl;
		return;
	}

	TFile* outfile = new TFile(output_rootfilename,"RECREATE");
	HistoManager *histoman = new HistoManager(outfile);
	histoman->loadHistoConfig("./HMConfig/_2body.HMconfig");
	DeterminePolygons(histoman);

	TMassTable fMassTable;
	//fMassTable.Init("/mnt/e/kinematics/IMMMA_Tool/threebody/masstable.dat");//uncomment this when on surface laptop -> eventually update this to be in a total config file
	//fMassTable.Init("/home/zmpur/IMMMA_Tool/threebody/masstable.dat");//uncomment this when on DESKTOP
	fMassTable.Init("../../config/masstable.dat");

	// //reaction:
	// std::vector<Reaction> reactions;
	// Reaction r;
	// r.beam.SetAll(3,"He",fMassTable.GetMassMeV("He",3));
	// r.target.SetAll(7,"Li",fMassTable.GetMassMeV("Li",7));
	// r.ejectile.SetAll(4,"He",fMassTable.GetMassMeV("He",4));
	// r.recoil.SetAll(6,"Li",fMassTable.GetMassMeV("Li",6));
	// r.breakup1.SetAll(4,"He",fMassTable.GetMassMeV("He",4));
	// r.breakup2.SetAll(2,"H",fMassTable.GetMassMeV("H",2));
	// r.beamEnergy = 7.5;
	// r.recoilExE = 1e-9;
	// reactions.push_back(r);

	// CMCMinMax Vcm1, Vcm2, KEcm1, KEcm2, ThetaCM1, ThetaCM2, PhiCM1, PhiCM2, Ecm, ThetaCMSum, PhiCMSep;

	// std::vector<IMMMA_Tool_3*> tools = prepareIMMMA_Tool_3s(reactions,true);
	// cout << "VCM1 (4He)  = " << tools[0]->GetExpectedVcm_bu1() << endl;
	// cout << "VCM2 (2H)   = " << tools[0]->GetExpectedVcm_bu2() << endl;
	// cout << "KECM1 (4He) = " << tools[0]->GetExpectedKEcm_bu1() << endl;
	// cout << "KECM2 (2H)  = " << tools[0]->GetExpectedKEcm_bu2() << endl;
	// cout << "ECM         = " << tools[0]->GetExpectedEcm() << endl;
	// cout << endl;

	//CUTS HERE WHEN APPLICABLE:


	TTree* kin2 = new TTree(ntpname, "2-body simulation tree");

	PHYSDATA physdata1, physdata2;//, physdata3, physdata4;		//holds the physics data for ejectile and recoil
	SABREDATA sabredata1, sabredata2;							//holds the sabre detection data for ejectile and recoil if detected by SABRE

	kin2->Branch("physdata1",&physdata1, "e/F:theta/F:phi/F");//branch to hold physics data for 
	kin2->Branch("physdata2",&physdata2, "e/F:theta/F:phi/F");

	kin2->Branch("sabredata1",&sabredata1,"detectorIndex/I:ring/I:wedge/I:theta/F:phi/F:ringEnergy/F:wedgeEnergy/F");
	kin2->Branch("sabredata2",&sabredata2,"detectorIndex/I:ring/I:wedge/I:theta/F:phi/F:ringEnergy/F:wedgeEnergy/F");

	//kin2->Branch();

	string line;
	//map<int,pair<int,int>> detector_hits;//ex: {404, {3,6}} means particle 4 in detector 4 in (r,w) = (3,6).
	//map<int, DetectorHit> detector_hits;//ex: detector_hits[402] = {3,6,4.56} means that particle 4 (first digit) was detected in SABRE2 (last digit) with energy 4.56 in (ring=3,wedge=6)
	bool hit1 = false, hit2 = false;

	cout << "Processing " << input_filename << "..." << endl;
	int count = 0;

	std::vector<std::string> eventLines;
	while(std::getline(infile,line)){
		if(line == "-1"){
			//cout << "got here, count = " << count << endl;	
			//event parsing
			PHYSDATA pd1, pd2;
			SABREDATA sd1, sd2;
			sd1.detectorIndex = -666;
			sd2.detectorIndex = -666;

			if(eventLines.size()==1){
				//only kinematics line
				parsePhysData(eventLines[0],pd1,pd2);
				fillKinHistos(histoman,pd1,pd2);
				//cout << "test" << endl;
			} else if(eventLines.size()==2){
				//kinematics1, particle 1 or particle 2
				parsePhysData(eventLines[0],pd1,pd2);
				parseSABREData(eventLines[1],sd1);
				sd1.theta = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].first;//wow this is ugly but it works
				sd1.phi = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].first;//likewise^
				fillKinHistos(histoman,pd1,pd2);
				fillSABREHistos(histoman,sd1,pd1);
				//cout << Form("%d\t%d\t%d\t%f\t%f\t%f\t%f",sd1.detectorIndex,sd1.ring,sd1.wedge,sd1.theta,sd1.phi,sd1.wedgeEnergy,sd1.ringEnergy) << endl;
				//cout << "test" << endl;
				Double_t exe = calculateSPS_ExE(pd1.e,pd1.theta,pd1.phi,fMassTable);
				histoman->getHisto2D("hSABRE_SabreRingEVsLi6ExE")->Fill(exe,sd1.ringEnergy);
				histoman->getHisto1D("hSABRE_SabreRingE")->Fill(sd1.ringEnergy);
				histoman->getHisto1D("hSPS_ExE")->Fill(exe);
				if(sd1.ring == 7){ }
				if(sd1.ring == 8){ }
				if(sd1.ring == 9){ }
			} else if(eventLines.size()==3){
				//kinematics1, particle 1 and particle 2
				parsePhysData(eventLines[0],pd1,pd2);
				parseSABREData(eventLines[1],sd1);
				sd1.theta = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].first;//wow this is ugly but it works
				sd1.phi = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].second;//likewise^
				parseSABREData(eventLines[2],sd2);
				sd2.theta = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first, sd2.wedge+offsets[sd2.detectorIndex].second}].first;//wow this is ugly but it works
				sd2.phi = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first, sd2.wedge+offsets[sd2.detectorIndex].second}].second;//likewise^
				fillKinHistos(histoman,pd1,pd2);
				fillSABREHistos(histoman,sd1,pd1);
				fillSABREHistos(histoman,sd2,pd2);
				Double_t exe = calculateSPS_ExE(pd1.e,pd1.theta,pd1.phi,fMassTable);
				histoman->getHisto2D("hSABRE_SabreRingEVsLi6ExE")->Fill(exe,sd1.ringEnergy);
				histoman->getHisto2D("hSABRE_SabreRingEVsLi6ExE")->Fill(exe,sd2.ringEnergy);
				histoman->getHisto1D("hSABRE_SabreRingE")->Fill(sd1.ringEnergy);
				histoman->getHisto1D("hSABRE_SabreRingE")->Fill(sd1.ringEnergy);
				histoman->getHisto1D("hSPS_ExE")->Fill(exe);

				//cout << "test" << endl;
			} else {
				cerr << "Warning: Unexpected eventLines.size() = " << eventLines.size() << endl;
			}

			sabredata1 = sd1;
			sabredata2 = sd2;
			physdata1 = pd1;
			physdata2 = pd2;

			//kin2->Write();

			eventLines.clear();
			count += 1;
			if(count%100000 == 0) cout << "Processed " << count << " events..." << endl;
		} else {
			eventLines.push_back(line);
		}
	}

	outfile->cd();
	kin2->Write();

	/*Update Histogram Axes Here*/
	UpdateHistoAxes(histoman);

	histoman->WriteAll(true);
	cout << endl;
	cout << "Processed " << count << " events." << endl;
	cout << "ROOT file saved to " << output_rootfilename << endl << endl;
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

	vector<TString> histonames = {
		"hSABRE0_ringEVSkinE",
		"hSABRE1_ringEVSkinE",
		"hSABRE2_ringEVSkinE",
		"hSABRE3_ringEVSkinE",
		"hSABRE4_ringEVSkinE",
		"hSABRE0_wedgeEVSkinE",
		"hSABRE1_wedgeEVSkinE",
		"hSABRE2_wedgeEVSkinE",
		"hSABRE3_wedgeEVSkinE",
		"hSABRE4_wedgeEVSkinE",
		"hSABRE0_ringThetaVSkinTheta",
		"hSABRE1_ringThetaVSkinTheta",
		"hSABRE2_ringThetaVSkinTheta",
		"hSABRE3_ringThetaVSkinTheta",
		"hSABRE4_ringThetaVSkinTheta",
		"hSABRE0_wedgePhiVSkinPhi",
		"hSABRE1_wedgePhiVSkinPhi",
		"hSABRE2_wedgePhiVSkinPhi",
		"hSABRE3_wedgePhiVSkinPhi",
		"hSABRE4_wedgePhiVSkinPhi"
	};

	for(const auto& name : histonames){
		histoman->getHisto2D(name)->GetXaxis()->SetTitle("SABRESim");
		histoman->getHisto2D(name)->GetYaxis()->SetTitle("Kin2mc");
		histoman->getHisto2D(name)->GetXaxis()->CenterTitle();
		histoman->getHisto2D(name)->GetYaxis()->CenterTitle();

		//histoman->getHisto2D(name)->SetOption("SQUARE");
	}

	histonames = {
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

#ifndef PLOT2MC_H
#define PLOT2MC_H

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "HistoManager.h"
#include "MassTable.h"
#include "TRandom3.h"
#include "Vec3.h"

struct PHYSDATA {double e, theta, phi;};
struct SABREDATA {int detectorIndex=-666, particleIndex=-666; double theta, phi, ringEnergy, wedgeEnergy, localx, localy; int ring, wedge;};

struct dp3Nucleus{
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

static const int numwedges = 8;
static const int numrings = 16;

static const float DEGRAD=0.017453293;

static const std::pair<int,int> offsets[] = {
	{112,40},	//detector0 {ringOffset,wedgeOffset}
	{96,32},	//detector1 {ringOffset,wedgeOffset}
	{80,16},	//detector2 {ringOffset,wedgeOffset}
	{64,24},	//detector3 {ringOffset,wedgeOffset}
	{48,0}		//detector4 {ringOffset,wedgeOffset}
};	

class plot2mc{
public:
	plot2mc(const std::string& outname);
	~plot2mc();

	void ProcessTXTOutput(std::vector<std::string> eventLines);
	void ProcessTXTOutput(const std::string& outputLines);

	void FillBeamSpotHisto(Vec3& reactionOrigin);

	void SaveAndWrite();

private:
	void FillKinematicsHistos(PHYSDATA& pd1, PHYSDATA& pd2);

	void FillSABREHistos(SABREDATA& sd1, PHYSDATA& pd1);

	bool ParsePhysData(const std::string& line, PHYSDATA& pd1, PHYSDATA& pd2);

	bool ParseSABREData(const std::string& line, SABREDATA& sd);

	double getSPSEnergy(double kinEnMeV, double sigmaMeV=0.015);

	double calculateSPS_ExE(double spsE, double spsTheta, double spsPhi);

	std::vector<std::string> splitLines(const std::string& ss_str){
		
		std::istringstream iss(ss_str);
		
		std::string line;

		std::vector<std::string> lines;

		while(std::getline(iss,line)){
			//if(!line.empty()) lines.push_back(line);
			lines.push_back(line);
		}

		return lines;
	}

	void readSingleAngleMap(std::ifstream& infile, std::map<std::pair<int,int>,std::pair<double,double>>& map);
	std::map<std::pair<int,int>,std::pair<double,double>> readAngleMaps();

	void UpdateHistoAxes();

	std::map<std::pair<int,int>, std::pair<double,double>> sabre_thetaphimap;

	TFile *histofile;
	HistoManager *histoman;
	MassTable *masstable;

	TRandom3 *rand_gen;
};


#endif//PLOT2MC_H
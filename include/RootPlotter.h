#ifndef ROOTPLOTTER_H
#define ROOTPLOTTER_H

#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "HistoManager.h"

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
		SetSym(newSym);
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

class RootPlotter{
public:
	RootPlotter(const std::string& outname);
	~RootPlotter();

	bool Init(const TString& configPath);

	void FillKinematics(int kinPID, double e, double theta, double phi);

	void FillSABRE(SABREDATA& sd, PHYSDATA& pd);

	void Analyze2Body(const char* input_filename, const char* output_rootfilename);

	void Analyze3


private:
	static const std::pair<int, int> offsets[];

	static const int numwedges = 8;
	static const int numrings = 16;

	HistoManager *histoman;



	std::pair<double,double>

};

#endif//ROOTPLOTTER_H
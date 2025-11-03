#ifndef ROOTPLOTTER_H
#define ROOTPLOTTER_H

#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "HistoManager.h"

struct PHYSDATA {double e, theta, phi;};
struct SABREDATA {int detectorIndex=-666, particleIndex=-666; double theta, phi, ringEnergy, wedgeEnergy, localx, localy; int ring, wedge;};

class RootPlotter{
public:
	RootPlotter(const std::string& outname);
	~RootPlotter();

	bool Init(const TString& configPath);

	void FillKinematics(int kinPID, double e, double theta, double phi);

	void FillSABRE(SABREDATA& sd, PHYSDATA& pd);




private:
	HistoManager *histoman;

};

#endif//ROOTPLOTTER_H
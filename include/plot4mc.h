#ifndef PLOT4MC_H
#define PLOT4MC_H

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
#include "structs.h"

class plot4mc {
public:
	plot4mc(const std::string& outname);
	~plot4mc();

	void ProcessTXTOutput(std::vector<std::string> eventLines);
	void ProcessTXTOutput(const std::string& outputLines);

	void FillBeamSpotHisto(Vec3& reactionOrigin);

	void FillStraggleHistos(double oldTheta, double oldPhi, double newTheta, double newPhi, double dTheta, double dPhi);

	void Fill_IMM(const CaseResult4& cr, double recoilGroundMassMeV);
	void Fill_MMM(const CaseResult4& cr, double recoilGroundMassMeV);

	bool FillTH1D(const TString& histoname, double value);
	bool FillTH2D(const TString& histoname, double valuex, double valuey);
	bool FillTH3D(const TString& histoname, double valuex, double valuey, double valuez);

	void SaveAndWrite();

private:
	void FillKinematicsHistos(PHYSDATA& pd1, PHYSDATA& pd2, PHYSDATA& pd3, PHYSDATA& pd4);

	void FillSABREHistos(SABREDATA& sd1, PHYSDATA& pd1);

	bool ParsePhysData(const std::string& line, PHYSDATA& pd1, PHYSDATA& pd2, PHYSDATA& pd3, PHYSDATA& pd4);

	bool ParseSABREData(const std::string& line, SABREDATA& sd);

	double getSPSEnergy(double kinEnMeV, double sigmaMeV=0.0015);

	double calculateSPS_ExE(double spsE, double spsTheta, double spsPhi);

	std::vector<std::string> splitLines(const std::string& ss_str){

		std::istringstream iss(ss_str);

		std::string line;

		std::vector<std::string> lines;

		while(std::getline(iss,line)){
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

#endif//PLOT4MC_H
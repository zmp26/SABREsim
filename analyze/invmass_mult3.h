#ifndef INVMASSMULT3_H
#define INVMASSMULT3_H

#include <map>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TDirectory.h"
#include "TLorentzVector.h"

class InvMass_Mult3{
private:
	std::map<TString, std::map<TString,TH1D*>> hMap;
	std::vector<TString> caseNames;

	double masses[3];
	double daughterMass, recoilMass;

	struct Perm { int i, j, k; };
	std::map<TString, Perm> pMap;

	const double DEGRAD = M_PI / 180.;
	const double RADDEG = 180. / M_PI;

	TFile *outfile;

public:
	InvMass_Mult3();
	~InvMass_Mult3();

	void Init(const char* output_filename);
	void SetMasses(double mass_bu1, double mass_bu2, double mass_bu3, double mass_recoil, double mass_daughter);

	void AnalyzeEvent(double E[3], double theta[3], double phi[3]);

	void CloseAndWrite();

};


#endif//INVMASSMULT3_H
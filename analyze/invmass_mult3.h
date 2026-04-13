#ifndef INVMASSMULT3_H
#define INVMASSMULT3_H

#include <map>
#include <vector>
#include <array>
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

	//define a "results" struct here:
	struct Results {
		double daughterIM, daughterExE, reconExE;

		double daughtervcm, daughterkecm, daughterthetacm, daughterphicm;
		double bu1vcm, bu1kecm, bu1thetacm, bu1phicm;
		double bu2vcm, bu2kecm, bu2thetacm, bu2phicm;
		double bu3vcm, bu3kecm, bu3thetacm, bu3phicm;

		double ecm1, ecm2;

		void Reset(){
			daughterIM = -666.;
			daughterExE = -666.;
			reconExE = -666.;

			daughtervcm = -666.;
			daughterkecm = -666.;
			daughterthetacm = -666.;
			daughterphicm = -666.;

			bu1vcm = -666.;
			bu1kecm = -666.;
			bu1thetacm = -666.;
			bu1phicm = -666.;

			bu2vcm = -666.;
			bu2kecm = -666.;
			bu2thetacm = -666.;
			bu2phicm = -666.;

			bu3vcm = -666.;
			bu3kecm = -666.;
			bu3thetacm = -666.;
			bu3phicm = -666.;

			ecm1 = -666.;
			ecm2 = -666.;
		};
	};

	Results caseResults[6];

	void ClearEventResults();

public:
	InvMass_Mult3();
	~InvMass_Mult3();

	void Init(const char* output_filename);
	void SetMasses(double mass_bu1, double mass_bu2, double mass_bu3, double mass_recoil, double mass_daughter);

	std::array<double,6> AnalyzeEvent(double E[3], double theta[3], double phi[3]);
	void FillEventHistograms();

	void CloseAndWrite();

};


#endif//INVMASSMULT3_H
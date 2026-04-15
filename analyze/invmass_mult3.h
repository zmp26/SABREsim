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
	std::vector<TString> permNames;

	double masses[3];
	double intermediateMass, recoilMass;

	struct Perm { int i, j, k; };
	std::map<TString, Perm> pMap;

	const double DEGRAD = M_PI / 180.;
	const double RADDEG = 180. / M_PI;

	TFile *outfile;

	//define a "results" struct here:
	struct Results {
		double intermediateIM, intermediateEx, reconEx;

		double intermediatevcm, intermediatekecm, intermediatethetacm, intermediatephicm;
		double frag1vcm, frag1kecm, frag1thetacm, frag1phicm;
		double frag2vcm, frag2kecm, frag2thetacm, frag2phicm;
		double frag3vcm, frag3kecm, frag3thetacm, frag3phicm;

		double ecm1, ecm2;

		void Reset(){
			intermediateIM = -666.;
			intermediateEx = -666.;
			reconEx = -666.;

			intermediatevcm = -666.;
			intermediatekecm = -666.;
			intermediatethetacm = -666.;
			intermediatephicm = -666.;

			frag1vcm = -666.;
			frag1kecm = -666.;
			frag1thetacm = -666.;
			frag1phicm = -666.;

			frag2vcm = -666.;
			frag2kecm = -666.;
			frag2thetacm = -666.;
			frag2phicm = -666.;

			frag3vcm = -666.;
			frag3kecm = -666.;
			frag3thetacm = -666.;
			frag3phicm = -666.;

			ecm1 = -666.;
			ecm2 = -666.;
		};
	};

	Results caseResults[6];

	void ClearEventResults();

	//expected CM constants:
	// struct ExpectedCM {
	// 	double Ecm1, Ecm2;
	// 	double vcm_frag1, vcm_frag2, vcm_frag3, vcm_intermediate;
	// 	double vkecm_frag1, vkecm_frag2, vkecm_frag3, vkecm_intermediate;
	// };
	// ExpectedCM expectedCMValues;
	// void SetExpectedCMValues();//called automatically at end of SetMasses()

public:
	InvMass_Mult3();
	~InvMass_Mult3();

	void Init(const char* output_filename);
	void SetMasses(double mass_frag1, double mass_frag2, double mass_frag3, double mass_recoil, double mass_intermediate);

	std::array<double,6> AnalyzeEvent(double E[3], double theta[3], double phi[3]);
	void FillEventHistograms();//fills all 6 cases together for the event
	void FillSelectCaseHistograms(int caseNum);//selectively fills a single case for the event

	void CloseAndWrite();

};


#endif//INVMASSMULT3_H
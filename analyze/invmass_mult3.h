#ifndef INVMASSMULT3_H
#define INVMASSMULT3_H

#include <map>
#include <vector>
#include <array>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "permHisto.h"

//Hypothesis4 struct to use in In
struct Hypothesis4 {
	std::string name;

	double mass_target;
	double mass_beam;
	double mass_ejectile;
	double mass_recoil;
	double mass_intermediate;

	double masses[3];//masses[0] = frag1, masses[1] = frag2, masses[2] = frag3

	double recoilEx;
	double intermediateEx;
	double intermediateExGate;
};

class InvMass_Mult3{
private:
	std::map<TString, std::map<TString,TH1*>> hMap;
	std::vector<TString> permNames;

	Hypothesis4 hypothesis;
	double masses[3];
	double intermediateMass, recoilMass, recoilEx;

	struct Perm { int i, j, k; };
	std::map<TString, Perm> pMap;

	const double DEGRAD = M_PI / 180.;
	const double RADDEG = 180. / M_PI;

	TFile *outfile;
	TTree *outtree;

	std::map<TString, permHisto*> groups_ungated;
	std::map<TString, permHisto*> groups_gated;

	double intermediateEx, intermediateExGate; //this holds the hypothesis of the intermediate/intermediate Ex and the gate (+/- due to width)

	//expected CM constants:
	struct ExpectedCM {
		double Ecm1, Ecm2;
		double vcm_frag1, vcm_frag2, vcm_frag3, vcm_intermediate;
		double kecm_frag1, kecm_frag2, kecm_frag3, kecm_intermediate;
	};

	//define a "results" struct here:
	struct Results {
		double intermediateIM, intermediateEx, reconEx;

		double intermediatevcm, intermediatekecm, intermediatethetacm, intermediatephicm;
		double frag1vcm, frag1kecm, frag1thetacm, frag1phicm;
		double frag2vcm, frag2kecm, frag2thetacm, frag2phicm;
		double frag3vcm, frag3kecm, frag3thetacm, frag3phicm;

		double ecm1, ecm2;

		ExpectedCM expected;

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

			expected.Ecm1 = -666.;
			expected.Ecm2 = -666.;
			expected.vcm_intermediate = -666.;
			expected.vcm_frag1 = -666.;
			expected.vcm_frag2 = -666.;
			expected.vcm_frag3 = -666.;
			expected.kecm_intermediate = -666.;
			expected.kecm_frag1 = -666.;
			expected.kecm_frag2 = -666.;
			expected.kecm_frag3 = -666.;
		};
	};

	Results caseResults[6];

	void ClearEventResults();

	ExpectedCM expectedCMValues;
	void SetExpectedCMValues(bool verbose=false);//called automatically at end of SetHypothesis()

public:
	InvMass_Mult3();
	~InvMass_Mult3();

	void Init(const char* output_filename);
	//void SetMasses(double mass_frag1, double mass_frag2, double mass_frag3, double mass_recoil, double mass_intermediate);
	void SetHypothesis(const Hypothesis4& hypo);

	std::array<double,6> AnalyzeEvent(double E[3], double theta[3], double phi[3], bool updateIntermediateEx=false);
	void FillEventHistograms();//fills all 6 cases together for the event
	void FillSelectCaseHistograms(int caseNum);//selectively fills a single case for the event

	void CloseAndWrite();

	void SetRecoilEx(double Ex) { recoilEx = Ex; hypothesis.recoilEx = Ex; SetExpectedCMValues(); }
	void SetIntermediateEx(double Ex) { intermediateEx = Ex; hypothesis.intermediateEx = Ex; SetExpectedCMValues(); }
	void SetIntermediateExGate(double ExGate) { intermediateExGate = ExGate; hypothesis.intermediateExGate = ExGate; }
	double GetDaughterEx() { return intermediateEx; }
	double GetDaughterExGate() { return intermediateExGate; }

};


#endif//INVMASSMULT3_H
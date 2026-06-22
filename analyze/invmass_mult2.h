#ifndef INVMASSMULT2_H
#define INVMASSMULT2_H

#include <map>
#include <vector>
#include <array>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "permHisto_mult2.h"

struct Hypothesis3 {

	std::string name;

	double mass_target;
	double mass_beam;
	double mass_ejectile;
	double mass_recoil;//ground state mass of decaying nucleus

	double masses[2];//masses[0] = frag1, masses[1] = frag2

};

enum gateIndicesMult2{
	NOCHECK_M2,
	FIRSTVALID_M2 = NOCHECK_M2,
	FRAG1VCMCHECK_M2,
	FRAG2VCMCHECK_M2,
	LASTVALID_M2 = FRAG2VCMCHECK_M2
};

class InvMass_Mult2{
private:
	std::map<TString, std::map<TString, TH1*>> hMap;
	std::vector<TString> permNames;

	Hypothesis3 hypothesis;
	double masses[2];
	double recoilMass, recoilEx;

	struct Perm{ int i, j; };
	std::map<TString, Perm> pMap;

	const double DEGRAD = M_PI/180.;
	const double RADDEG = 180./M_PI;

	TFile *outfile;
	TTree *outtree;

	int gate1_index;
	int gate2_index;

	std::pair<double, double> gate1minmax;
	std::pair<double, double> gate2minmax;

	std::map<TString, permHisto_mult2*> groups_ungated;
	std::map<TString, permHisto_mult2*> groups_gated;

	struct ExpectedCM_M2{
		double Ecm;
		double vcm_frag1, kecm_frag1;
		double vcm_frag2, kecm_frag2;
	};

	struct Results {
		double recoilIM, reconEx;
		double frag1vcm, frag1kecm, frag1thetacm, frag1phicm;
		double frag1Comp[3];
		double frag2vcm, frag2kecm, frag2thetacm, frag2phicm;
		double frag2Comp[3];

		double ecm;
		double boost[3];
		double relLabAngle_frag1frag2;

		double exp_ecm;
		double exp_f1VCM, exp_f1KECM;
		double exp_f2VCM, exp_f2KECM;

		bool permPasses = false;
		TString permName;
		ExpectedCM_M2 expected;

		void Reset(){
			permName = "NONE";

			recoilIM = -666.;
			reconEx = -666.;

			frag1vcm = -666.;
			frag1kecm = -666.;
			frag1thetacm = -666.;
			frag1phicm = -666.;

			frag2vcm = -666.;
			frag2kecm = -666.;
			frag2thetacm = -666.;
			frag2phicm = -666.;

			ecm = -666.;
			relLabAngle_frag1frag2 = -666.;
			permPasses = false;

			for(int i=0; i<3; i++){
				boost[i] = -666.;
				frag1Comp[i] = -666.;
				frag2Comp[i] = -666.;
			}

			expected.Ecm = -666.;
			expected.vcm_frag1 = -666.;
			expected.vcm_frag2 = -666.;
			expected.kecm_frag1 = -666.;
			expected.kecm_frag2 = -666.;
		}
	};

	Results caseResults[2];

	void ClearEventResults();

	ExpectedCM_M2 expectedCMValues;
	void SetExpectedCMValues(bool verbose = false);

	TH1D* hPermCounter;
	TH1D* hPermCounter_gated;
	TH1D *hSortedPermutations;
	TH1D *hSortedIMRecEx;

public:
	InvMass_Mult2();
	~InvMass_Mult2();

	void Init(const char* output_filename);
	void SetHypothesis(const Hypothesis3& hypo);

	std::array<double, 2> AnalyzeEvent(double E[2], double theta[2], double phi[2]);
	void FillEventHistograms(double SPS_Ex);
	void FillSelectCaseHistograms(int caseNum, double SPS_Ex);

	void FillSortedHisto(double SPS_Ex);

	void FillGatedEventHistograms(double SPS_Ex);
	void FillSelectGatedCaseHistograms(int caseNum, double SPS_Ex);

	void FillTree();
	void FillPermCounter(bool gated=false);
	void CloseAndWrite();

	void SetRecoilEx(double Ex){ recoilEx = Ex; SetExpectedCMValues(); };

	void SetGate1(int index){ if(FIRSTVALID_M2 <= index && LASTVALID_M2 >= index) gate1_index = index; }
	void SetGate2(int index){ if(FIRSTVALID_M2 <= index && LASTVALID_M2 >= index) gate2_index = index; }
	void SetGate1MinMax(std::pair<double, double> minmax){ gate1minmax = minmax; };
	void SetGate2MinMax(std::pair<double,double> minmax){ gate2minmax = minmax; };

	bool CheckGate1(double val){return(val >= gate1minmax.first && val <= gate1minmax.second);}
	bool CheckGate2(double val){return(val >= gate2minmax.first && val <= gate2minmax.second);}

	int CountPermPasses();
};

#endif//INVMASSMULT2_H
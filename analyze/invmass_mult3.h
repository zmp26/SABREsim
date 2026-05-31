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
#include "permHisto_mult3.h"

//Hypothesis4 struct
struct Hypothesis4 {
	std::string name;

	double mass_target;
	double mass_beam;
	double mass_ejectile;
	double mass_recoil;
	double mass_intermediate;

	double masses[3];//masses[0] = frag1, masses[1] = frag2, masses[2] = frag3

	// double recoilEx;
	// double intermediateEx;
	// double intermediateExGate;
};

enum gateIndices {
	NOCHECK,
	FIRSTVALID = NOCHECK,
	INTEXCHECK,
	INTVCMCHECK,
	FRAG1VCMCHECK,
	LASTVALID = FRAG1VCMCHECK
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

	/*
		gate1/2_index mapping:
			0 	=	no gate check (just pass true to if statement)
			1	=	intermediateExCheck
			2	=	intermediateVcmCheck
			3 	=	frag1VcmCheck

		Default is gate1_index = 1, gate2_index = 0 for basic intermediateEx gate

		This should be set via invmass_mult3::SetGate1/2() before beginning analysis of events

		**SEE ENUM DECLARED OUTSIDE OF CLASS AT TOP OF FILE**
	*/

	int gate1_index;
	int gate2_index;
	int gate3_index;

	std::pair<double,double> gate1minmax;
	std::pair<double,double> gate2minmax;
	std::pair<double,double> gate3minmax;

	std::map<TString, permHisto_mult3*> groups_ungated;
	std::map<TString, permHisto_mult3*> groups_gated;

	//double intermediateEx, intermediateExGate; //this holds the hypothesis of the intermediate/intermediate Ex and the gate (+/- due to width)
	double intermediateEx;//, intermediateEmin, intermediateEmax;

	//expected CM constants:
	struct ExpectedCM {
		double Ecm1, Ecm2;
		// double vcm_frag1, vcm_frag2, vcm_frag3, vcm_intermediate;
		// double kecm_frag1, kecm_frag2, kecm_frag3, kecm_intermediate;
		double vcm_intermediate, kecm_intermediate;
		double vcm_frag1, kecm_frag1;
		double vcm_frag2, kecm_frag2;
		double vcm_frag3, kecm_frag3;
	};

	//define a "results" struct here:
	struct Results {
		TString permName;

		double intermediateIM, intermediateEx, reconEx;

		double intermediatevcm, intermediatekecm, intermediatethetacm, intermediatephicm;
		double intermediateComp[3];
		double frag1vcm, frag1kecm, frag1thetacm, frag1phicm;
		double frag1Comp[3];
		double frag2vcm, frag2kecm, frag2thetacm, frag2phicm;
		double frag2Comp[3];
		double frag3vcm, frag3kecm, frag3thetacm, frag3phicm;
		double frag3Comp[3];

		double ecm1, ecm2;

		double boost1[3], boost2[3];

		double relLabAngle_intfrag1, relLabAngle_frag2frag3;

		//double IM2_intfrag1, IM2_frag2frag3;
		//double IM2_01, IM2_12, IM2_20;
		double IM2_int;

		ExpectedCM expected;

		bool permPasses = false;

		void Reset(){
			permName = "NONE";

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

			relLabAngle_intfrag1 = -666.;
			relLabAngle_frag2frag3 = -666.;

			permPasses = false;

			for(int i=0; i<3; i++){
				boost1[i] = -666.;
				boost2[i] = -666.;
				intermediateComp[i] = -666.;
				frag1Comp[i] = -666.;
				frag2Comp[i] = -666.;
				frag3Comp[i] = -666.;
			}

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

	//"correct" permutation histograms:
	TH1D* hPermCounter;
	TH1D* hPermCounter_gated;

	//histogram(s) for caseResult entry with recoilEx nearest SPS provided value:
	TH2D* hSortedIntermediateExIMvsSPS;
	TH1D* hSortedPermutations;
	TH1D* hSortedIMRecEx;
	TH1D* hSortedIMRecEx_gate8Be;
	TH1D* hSortedIMRecEx_gate5Li;

public:
	InvMass_Mult3();
	~InvMass_Mult3();

	void Init(const char* output_filename);
	//void SetMasses(double mass_frag1, double mass_frag2, double mass_frag3, double mass_recoil, double mass_intermediate);
	void SetHypothesis(const Hypothesis4& hypo);

	std::array<double,6> AnalyzeEvent(double E[3], double theta[3], double phi[3], bool updateIntermediateEx=false);
	void FillEventHistograms(double SPS_Ex);//fills all 6 cases together for the event - ungated only
	void FillSelectCaseHistograms(int caseNum, double SPS_Ex);//selectively fills a single case for the event - ungated only

	void FillGatedEventHistograms(double SPS_Ex);//fills all 6 cases together for the event - gated only (assumes check done on user side!)
	void FillSelectGatedCaseHistograms(int caseNum, double SPS_Ex);//selectively fills a single case for the event - gated only (assumes check done on user side!)

	void FillTree();

	void FillPermCounter(bool gated=false);

	void FillSortedHisto(double SPS_Ex);

	void CloseAndWrite();

	void SetRecoilEx(double Ex) { recoilEx = Ex; SetExpectedCMValues(); }// hypothesis.recoilEx = Ex; SetExpectedCMValues(); }
	void SetIntermediateEx(double Ex) { intermediateEx = Ex; SetExpectedCMValues(); }// hypothesis.intermediateEx = Ex; SetExpectedCMValues(); }
	//void SetIntermediateExGate(double ExGate) { intermediateExGate = ExGate; }// hypothesis.intermediateExGate = ExGate; }
	// double GetIntermediateEx() { return intermediateEx; }
	// double GetIntermediateExGate() { return intermediateExGate; }
	void SetGate1(int index) { if(FIRSTVALID <= index && LASTVALID >= index) gate1_index = index; }
	void SetGate2(int index) { if(FIRSTVALID <= index && LASTVALID >= index) gate2_index = index; }
	void SetGate3(int index) { if(FIRSTVALID <= index && LASTVALID >= index) gate3_index = index; }
	void SetGate1MinMax(std::pair<double,double> minmax) { gate1minmax = minmax; }
	void SetGate2MinMax(std::pair<double,double> minmax) { gate2minmax = minmax; }
	void SetGate3MinMax(std::pair<double,double> minmax) { gate3minmax = minmax; }

	bool CheckGate1(double val) { return(val >= gate1minmax.first && val <= gate1minmax.second); }
	bool CheckGate2(double val) { return(val >= gate2minmax.first && val <= gate2minmax.second); }
	bool CheckGate3(double val) { return(val >= gate3minmax.first && val <= gate3minmax.second); }

	int CountPermPasses();

};


#endif//INVMASSMULT3_H
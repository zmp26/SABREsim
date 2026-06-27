#ifndef MISSMASSMULT2_H
#define MISSMASSMULT2_H

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
#include "permHistoMM_mult2.h"

/*
	missmass_mult2::Hypothesis4MM differs from invmass_mult3::Hypothesis4
	only by the addition of the beamEnergyMeV variable. This is necessary
	for missmass_mult2 case because the missing particle is reconstructed
	using the kinematical infomratoin from the assumed reaction (beam,
	target, ejectile).
*/
struct Hypothesis4MM {
	std::string name;

	double mass_target;
	double mass_beam;
	double mass_ejectile;
	double mass_recoil;
	double mass_intermediate;

	double beamEnergyMeV;

	double masses[3];//masses[0]=frag1, masses[1]=frag2, masses[2]=frag3
};

enum gateIndices {
	NOCHECK,
	FIRSTVALID = NOCHECK,
	INTEXCHECK,
	INTVCMCHECK,
	FRAG1VCMCHECK,
	LASTVALID = FRAG1VCMCHECK
};

class MissMass_Mult2{
private:
	std::map<TString, std::map<TString, TH1*>> hMap;
	std::vector<TString> permNames;

	Hypothesis4MM hypothesis;
	double masses[3];
	double beamMass, targetMass, ejectileMass;
	double beamEnergyMeV;
	double intermediateMass, recoilMass, recoilEx;

	struct Perm { int i, j, k; };
	std::map<TString, Perm> pMap;

	const double DEGRAD = M_PI/180.;
	const double RADDEG = 180./M_PI;
	const double nucleonMassMeV = (938.27 + 939.57)/2.;

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

	int gate1_index, gate2_index;

	std::pair<double,double> gate1minmax, gate2minmax;

	std::map<TString, permHistoMM_mult2*> groups_ungated;
	std::map<TString, permHistoMM_mult2*> groups_gated;

	double intermediateEx;

	struct ExpectedCM{
		double Ecm1, Ecm2;
		double vcm_intermediate, kecm_intermediate;
		double vcm_frag1, kecm_frag1;
		double vcm_frag2, kecm_frag2;
		double vcm_frag3, kecm_frag3;
	};

	struct ResultsMM{

		double intermediateIM, intermediateEx, reconEx;
		double intermediatevcm, intermediatekecm, intermediatethetacm, intermediatephicm;
		double intermediateComp[3];
		double frag1vcm, frag1kecm, frag1thetacm, frag1phicm;
		double frag1Comp[3];
		double frag2vcm, frag2kecm, frag2thetacm, frag2phicm;
		double frag2Comp[3];
		double frag3vcm, frag3kecm, frag3thetacm, frag3phicm;
		double frag3Comp[3];

		double missingmass;
		double SABREsumE;

		double ecm1, ecm2;
		double boost1[3], boost2[3];
		double relLabAngle_intfrag1, relLabAngle_frag2frag3;
		double IM2_int;

		double exp_ecm1, exp_ecm2;
		double exp_imVCM, exp_imKECM;
		double exp_f1VCM, exp_f1KECM;
		double exp_f2VCM, exp_f2KECM;
		double exp_f3VCM, exp_f3KECM;

		double catania_x, catania_y;

		bool permPasses = false;
		
		TString permName;
		ExpectedCM expected;

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

			missingmass = -666.;
			SABREsumE = -666.;

			ecm1 = -666.;
			ecm2 = -666.;

			relLabAngle_intfrag1 = -666.;
			relLabAngle_frag2frag3 = -666.;

			catania_x = -666.;
			catania_y = -666.;

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

	ResultsMM caseResults[6];

	void ClearEventResults();

	ExpectedCM expectedCMValues;
	void SetExpectedCMValues(bool verbose=false);

	TH1D* hPermCounter;
	TH1D* hPermCounter_gated;

	//histogram(s) for caseResult entry with recoilEx nearest SPS provided value:
	TH2D* hSortedIntermediateExIMvsSPS;
	TH2D* hSortedCataniaPlot;
	TH1D* hSortedPermutations;
	TH1D* hSortedIMRecEx;
	TH1D* hSortedIMRecEx_gate8Be;
	TH1D* hSortedIMRecEx_gate5Li;
	TH2D* hSABRESumE_vs_ExSPS;

public:
	MissMass_Mult2();
	~MissMass_Mult2();

	void Init(const char* output_filename);
	void SetHypothesis(const Hypothesis4MM& hypo);

	std::array<double,6> AnalyzeEvent(double E[2], double theta[2], double phi[2], double SPSE, double SPSTheta, double SPSPhi, bool updateIntermediateEx=true);
	void FillEventHistograms(double SPS_Ex);
	void FillSelectCaseHistograms(int caseNum, double SPS_Ex);

	void FillGatedEventHistograms(double SPS_Ex);
	void FillSelectGatedCaseHistograms(int caseNum, double SPS_Ex);

	void FillTree();

	void FillPermCounter(bool gated=false);

	void FillSortedHisto(double SPS_Ex);

	void FillSABREvsSPSHisto(double SPS_Ex, double SABREsumE);

	void CloseAndWrite();

	void SetRecoilEx(double Ex) { recoilEx = Ex; SetExpectedCMValues(); }
	void SetIntermediateEx(double Ex) { intermediateEx = Ex; SetExpectedCMValues(); }

	void SetGate1(int index) { if(FIRSTVALID <= index && LASTVALID >= index) gate1_index = index; }
	void SetGate2(int index) { if(FIRSTVALID <= index && LASTVALID >= index) gate2_index = index; }

	void SetGate1MinMax(std::pair<double,double> minmax) { gate1minmax = minmax; }
	void SetGate2MinMax(std::pair<double,double> minmax) { gate2minmax = minmax; }

	bool CheckGate1(double val) { return(val >= gate1minmax.first && val <= gate1minmax.second); }
	bool CheckGate2(double val) { return(val >= gate2minmax.first && val <= gate2minmax.second); }

	int CountPermPasses();

};

#endif//MISSMASSMULT2_H
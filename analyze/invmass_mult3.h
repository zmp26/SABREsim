#ifndef INVMASSMULT3_H
#define INVMASSMULT3_H

#include <map>
#include <vector>
#include <array>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "permHisto_mult3.h"
#include "CutHandler.h"

/*
	invmass_mult3::Hypothesis4 differs from missmass_mult2::Hypothesis4MM
	only by the removal of the beamEnergyMeV variable. This is necessary
	for missmass_mult2 case because the missing particle is reconstructed
	using the kinematical infomratoin from the assumed reaction (beam,
	target, ejectile). However, it is not strictly needed in the invariant
	mass case (N case).
*/
struct Hypothesis4 {
	std::string name;

	double mass_target;
	double mass_beam;
	double mass_ejectile;
	double mass_recoil;
	double mass_intermediate;

	double beamEnergyMeV;

	//note that kin4mc treats masses[0] as the first decay step and masses[2] as the second decay step, leaving masses[3] the "daughter"
	//For example, 5Li -> p + a means that 5Li is proton decaying and 5Li -> a + p is 5Li undergoing alpha decay
	//Yes, these are energetically identical in the above case
	//However, one can imagine heaver resonances resulting in a heavier, but still stable, resonance decay particles for masses[3]
	//For this reason, we will calculate an "above masses[2] threshold" to generate a plot of intermediate energy above thresohld
	//(This is essentially just another shifted histogram of the true invariant mass spectrum)
	//We will do the same for the initial resonance (the recoil from reaction that decays into, for example, the 5Li above)
	double masses[3];//masses[0] = frag1, masses[1] = frag2, masses[2] = frag3

	// double recoilEx;
	// double intermediateEx;
	// double intermediateExGate;
};

class InvMass_Mult3{
private:
	std::map<TString, std::map<TString,TH1*>> hMap;
	std::vector<TString> permNames;

	Hypothesis4 hypothesis;
	double masses[3];
	double intermediateMass, recoilMass, recoilEx;
	double beamEnergyMeV;

	struct Perm { int i, j, k; };
	std::map<TString, Perm> pMap;

	const double DEGRAD = M_PI / 180.;
	const double RADDEG = 180. / M_PI;

	TFile *outfile;
	TTree *outtree;
	bool writeTree = false;

	std::map<TString, permHisto_mult3*> groups_ungated;
	std::map<TString, permHisto_mult3*> groups_gated;

	//double intermediateEx, intermediateExGate; //this holds the hypothesis of the intermediate/intermediate Ex and the gate (+/- due to width)
	double intermediateEx;//, intermediateEmin, intermediateEmax;

	double intermediateM1Threshold;
	double recoilM0Threshold;

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

		double SPSEnergy, SPSTheta, SPSPhi, SPS_Ex;

		double intermediateIM, intermediateEx, reconEx;
		double intermediateEnergyAboveM1Thresh;
		double recoilEnergyAboveM0Thresh;
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
		double relLabAngle_intfrag1, relLabAngle_frag1frag2, relLabAngle_frag2frag3;
		double IM2_int;

		double Theta2h;
		double CosTheta2h;

		double exp_ecm1, exp_ecm2;
		double exp_imVCM, exp_imKECM;
		double exp_f1VCM, exp_f1KECM;
		double exp_f2VCM, exp_f2KECM;
		double exp_f3VCM, exp_f3KECM;

		double m12sq, m23sq;

		double PLabTotal_ResDecayParticles;
		double PLabTotal_Beam;
		double PLabTotal_Ejectile;
		double MissingMomentumComp[3];
		double MissingMomentumMag;
		double E1meas, E2meas, E3meas;
		double Theta1meas, Theta2meas, Theta3meas;

		bool permPasses = false;
		
		TString permName;
		ExpectedCM expected;

		void Reset(){
			permName = "NONE";

			SPSEnergy = -666.;
			SPSTheta = -666.;
			SPSPhi = -666.;
			SPS_Ex = -666.;

			intermediateIM = -666.;
			intermediateEx = -666.;
			reconEx = -666.;

			intermediateEnergyAboveM1Thresh = -666.;
			recoilEnergyAboveM0Thresh = -666.;

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
			relLabAngle_frag1frag2 = -666.;
			relLabAngle_frag2frag3 = -666.;
			IM2_int = -666.;

			Theta2h = -666.;
			CosTheta2h = -666.;

			m12sq = -666.;
			m23sq = -666.;

			PLabTotal_Beam = -666.;
			PLabTotal_ResDecayParticles = -666.;
			PLabTotal_Ejectile = -666.;
			E1meas = -666.;
			E2meas = -666.;
			E3meas = -666.;
			Theta1meas = -666.;
			Theta2meas = -666.;
			Theta3meas = -666.;
			MissingMomentumMag = -666.;

			permPasses = false;

			for(int i=0; i<3; i++){
				boost1[i] = -666.;
				boost2[i] = -666.;
				intermediateComp[i] = -666.;
				frag1Comp[i] = -666.;
				frag2Comp[i] = -666.;
				frag3Comp[i] = -666.;
				MissingMomentumComp[i] = -666.;
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
	TH2D* hSortedDalitz;
	TH1D* hSortedPermutations;
	TH1D* hSortedIMRecEx;
	TH1D* hSortedIMRecEx_gate8Be;
	TH1D* hSortedIMRecEx_gate5Li;
	TH2D* hSABRESumE_vs_ExSPS;

	CutHandler fCutHandler;

	std::map<TString, double> PackMetrics1D(const Results& res) const;
	std::map<TString, std::pair<double, double>> PackPoints2D(const Results& res) const;

public:
	InvMass_Mult3();
	~InvMass_Mult3();

	void Init(const char* output_filename);
	//void SetMasses(double mass_frag1, double mass_frag2, double mass_frag3, double mass_recoil, double mass_intermediate);
	void SetHypothesis(const Hypothesis4& hypo);

	std::array<double,6> AnalyzeEvent(double E[3], double theta[3], double phi[3], double SPS_E, double SPSTheta, double SPSPhi, double SPS_Ex, bool updateIntermediateEx=false);
	void FillEventHistograms(double SPS_Ex);//fills all 6 cases together for the event - ungated only

	void FillGatedEventHistograms(double SPS_Ex);//fills all 6 cases together for the event - gated only (assumes check done on user side!)

	void FillTree();

	void FillPermCounter(bool gated=false);

	void FillSortedHisto(double SPS_Ex);

	void FillSABRESumEVsSPS_Ex(double SPS_Ex, double SABREsumE);

	void SetWriteTree(bool val){ writeTree = val; }
	bool GetWriteTree() { return writeTree; }

	void CloseAndWrite();

	void SetRecoilEx(double Ex) { recoilEx = Ex; SetExpectedCMValues(); }// hypothesis.recoilEx = Ex; SetExpectedCMValues(); }
	void SetIntermediateEx(double Ex) { intermediateEx = Ex; SetExpectedCMValues(); }// hypothesis.intermediateEx = Ex; SetExpectedCMValues(); }

	CutHandler& GetCutHandler() { return fCutHandler; }

	int CountPermPasses();
};


#endif//INVMASSMULT3_H
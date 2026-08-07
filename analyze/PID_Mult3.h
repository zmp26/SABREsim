#ifndef PID_MULT3_H
#define PID_MULT3_H

#include <array>
#include <vector>
#include <string>
#include <cmath>

#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"

// Streamlined hypothesis: Only global reaction and final state particle masses
struct PIDHypothesis {
	std::string name;

	// Reaction components (in MeV or MeV/c^2)
	double mass_target{0.0};     
	double mass_beam{0.0};       
	double mass_ejectile{0.0};   
	double beamEnergyMeV{0.0};

	// Final state particles: [0], [1], [2]
	double final_masses[3]{0.0, 0.0, 0.0};//holds the masses of final state particles [0], [1], [2] (for calculations)
	TString final_particles[3]{"","",""};//holds the symbol of final state particles [0], [1], [2] (for histogram labeling)
	//note there is not check for consistency between final_particles and final_masses, it assumes it is correct! You have been warned!
};

struct PIDResult { 
	// Index mapping: hit_indices[i] gives the detector hit index assigned to final_masses[i]
	std::array<int, 3> hit_indices = {-1, -1, -1};

	int bestChi2Index = -1;
	double bestChi2 = 1e9;

	double missing_px = 0.0;
	double missing_py = 0.0;
	double missing_pz = 0.0;
	double missing_E  = 0.0;
	double missing_Pmag = 0.0;

	std::array<double, 6> permChi2s;

	bool passesCut = false;

	void Reset() {
		hit_indices = {-1, -1, -1};
		bestChi2 = 1e9;
		missing_px = missing_py = missing_pz = missing_E = missing_Pmag = 0.0;
		bestChi2Index = -1;
		permChi2s.fill(1e9);
		passesCut = false;
	}
};

class PID_Mult3 {
private:
	PIDHypothesis fHypothesis;

	// Full set of 6 permutations for 3 hits (i, j, k)
	const std::array<std::array<int, 3>, 6> allPerms = {{
		{0, 1, 2}, {0, 2, 1},
		{1, 0, 2}, {1, 2, 0},
		{2, 0, 1}, {2, 1, 0}
	}};

	// Resolution parameters
	double fSigmaE = 0.050;		// MeV
	double fSigmaTheta = 1.0;	// degrees
	double fSigmaPhi = 1.0;		// degrees

	// SPS focal plane ejectile resolutions
	double fSPSSigmaE = 0.015;		// MeV (15 keV)
	double fSPSSigmaTheta = 1.0;	// degrees
	double fSPSSigmaPhi = 1.0;		// degrees

	// max chi2 to consider "success"
	double fMaxChi2Cut = 10.0;

	TDirectory* outdir{nullptr};
	TH1D* hBestChi2{nullptr};
	TH1D* hBestPermutation{nullptr};
	TH2D* hChi2_BestVsNext{nullptr};
	TH2D* h2Chi2ByPermutation{nullptr};
	TH1D* hProtonPicks{nullptr};//which hit index is picked as the proton

	TLorentzVector Build4Vector(double EMeV, double thetadeg, double phideg, double massMeV) const;
    
	std::array<std::array<double, 3>, 3> GetMomentumCovariance(
		double E, double theta_deg, double phi_deg, double mass, 
		double sigE, double sigTheta_deg, double sigPhi_deg) const;

	double ComputeCovarianceChi2(const std::array<std::array<double,3>,3>& Cov, const TLorentzVector& residual) const;

public:
	PID_Mult3();
	~PID_Mult3();

	void SetHypothesis(const PIDHypothesis& hypo);
	void SetResolution(double sigE, double sigTheta, double sigPhi);
	void SetSPSResolution(double sigE, double sigTheta, double sigPhi);
	void SetChi2Cut(double maxChi2) { fMaxChi2Cut = maxChi2; }

	// Optional ROOT logging
	void InitDiagnostics(TDirectory* targetDir);

	// Core PID solver
	PIDResult EvaluateEvent(const double E[3], const double theta[3], const double phi[3],
							double SPS_E, double SPSTheta, double SPSPhi);

	// Single-permutation Chi2 evaluator
	double ComputePermutationChi2(const std::array<int, 3>& perm,
								  const double E[3], const double theta[3], const double phi[3],
								  double SPS_E, double SPSTheta, double SPSPhi,
								  const TLorentzVector& P_recoil_expected,
								  TLorentzVector& P_residual_out) const;

	const PIDHypothesis& GetHypothesis() const { return fHypothesis; }
};

#endif // PID_MULT3_H
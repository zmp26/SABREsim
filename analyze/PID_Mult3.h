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

#include "PID_structs.h"

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
	TH2D* h2Chi2DifByPermutation{nullptr};
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
	PIDResult_Mult3 EvaluateEvent(const double E[3], const double theta[3], const double phi[3],
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
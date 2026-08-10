#ifndef PID_MULT2_H
#define PID_MULT2_H

#include <array>
#include <vector>
#include <string>
#include <cmath>

#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

#include "PID_structs.h"

class PID_Mult2 {
public:
	PID_Mult2();
	virtual ~PID_Mult2();

	void SetHypothesis(const PIDHypothesis& hypo);
	void SetResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg);
	void SetSPSResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg);
	void SetChi2Cut(double cut) { fMaxChi2Cut = cut; }

	void InitDiagnostics(TDirectory* targetDir);

	PIDResult_Mult2 EvaluateEvent(const double E[2], const double theta[2], const double phi[2],
								  const double SPS_E, const double SPSTheta, const double SPSPhi);

private:
	PIDHypothesis fHypothesis;

	double fSigmaE{0.};
	double fSigmaTheta{0.};
	double fSigmaPhi{0.};

	double fSPSSigmaE{0.};
	double fSPSSigmaTheta{0.};
	double fSPSSigmaPhi{0.};

	double fMaxChi2Cut{10.};

	TDirectory* outdir{nullptr};

	TH1D* hBestChi2{nullptr};
    TH1D* hBestPermutation{nullptr};
    TH1D* hMissingMassBest{nullptr};
    TH1D* hSigmaMBest{nullptr};
    TH2D* hChi2_BestVsNext{nullptr};
    TH2D* h2Chi2ByPermutation{nullptr};
    TH2D* h2Chi2DifByPermutation{nullptr};

    TLorentzVector Build4Vector(double EMeV, double thetadeg, double phideg, double massMeV) const;

    std::array<std::array<double,4>,4> GetP4Covariance(double EMeV, double thetadeg, double phideg, double massMeV,
															  double sigE, double sigThetaDeg, double sigPhiDeg) const;

    double ComputePermutationChi2(const std::array<int,3>& perm,
    							  const double E[2], const double theta[2], const double phi[2],
    							  double SPS_E, double SPSTheta, double SPSPhi,
    							  const TLorentzVector& P_recoil_expected,
    							  TLorentzVector& P_missing_out,
    							  double& sigma_M_out) const;

    const std::array<std::array<int,3>,6> allPerms = {{
    	{0,1,2}, {0,2,1},
    	{1,0,2}, {1,2,0},
    	{2,0,1}, {2,1,0}
    }};

};

#endif // PID_MULT2_H
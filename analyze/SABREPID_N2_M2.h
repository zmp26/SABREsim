#ifndef SABREPID_N2_M2_H
#define SABREPID_N2_M2_H

#include <array>
#include <cmath>

#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"

#include "PID_structs.h"

class SABREPID_N2_M2{
public:
	SABREPID_N2_M2();
	virtual ~SABREPID_N2_M2();

	void InitDiagnostics(TDirectory* targetDir);

	void SetHypothesis(const PIDHypothesis_N2& hypo);
	void SetResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg);
	void SetSPSResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg);
	void SetChi2Cut(double cut) { fMaxChi2Cut = cut; }

	PIDResult_N2_M2 EvaluateEvent(const double E[2], const double theta[2], const double phi[2], double SPS_E, double SPSTheta, double SPSPhi);

	const PIDHypothesis_N2& GetHypothesis() const { return fHypothesis; }

private:
	PIDHypothesis_N2 fHypothesis;

	double fSigmaE{0.05};
	double fSigmaTheta{1.};
	double fSigmaPhi{1.};

	double fSPSSigmaE{0.015};
	double fSPSSigmaTheta{1.};
	double fSPSSigmaPhi{1.};

	double fMaxChi2Cut{10.};

	const std::array<std::array<int,2>,2> allPerms = {{ {0,1}, {1,0} }};

	TDirectory* outdir{nullptr};

	TH1D* hBestChi2{nullptr};
	TH1D* hBestPermutation{nullptr};
	TH2D* hChi2_BestVsNext{nullptr};
	TH2D* h2Chi2ByPermutation{nullptr};
	TH2D* h2Chi2DifByPermutation{nullptr};

	TLorentzVector Build4Vector(double EMeV, double thetaDeg, double phiDeg, double massMeV) const;

	std::array<std::array<double,3>,3> GetMomentumCovariance(double E, double thetaDeg, double phiDeg, double mass, double sigE, double sigThetaDeg, double sigPhiDeg) const;

	double ComputeCovarianceChi2(const std::array<std::array<double,3>,3>& Cov, const TLorentzVector& residual) const;

	double ComputePermutationChi2(const std::array<int,2>& perm, const double E[2], const double theta[2], const double phi[2], double SPS_E, double SPSTheta, double SPSPhi, const TLorentzVector& P_recoil_expected, TLorentzVector& P_residual_out) const;
	
};

#endif //SABREPID_N2_M2_H
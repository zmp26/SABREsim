#include "SABREPID_N2_M2.h"
#include "PID_structs.h"

#include <cmath>

#include "TString.h"

SABREPID_N2_M2::SABREPID_N2_M2() = default;

SABREPID_N2_M2::~SABREPID_N2_M2() = default;

void SABREPID_N2_M2::SetHypothesis(const PIDHypothesis_N2& hypo){
	fHypothesis = hypo;
}

void SABREPID_N2_M2::SetResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg){
	fSigmaE = sigEMeV;
	fSigmaTheta = sigThetaDeg;
	fSigmaPhi = sigPhiDeg;
}

void SABREPID_N2_M2::SetSPSResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg){
	fSPSSigmaE = sigEMeV;
	fSPSSigmaTheta = sigThetaDeg;
	fSPSSigmaPhi = sigPhiDeg;
}




TLorentzVector SABREPID_N2_M2::Build4Vector(double EMeV, double thetadeg, double phideg, double massMeV) const {
	double p = std::sqrt(EMeV*EMeV + 2.*EMeV*massMeV);
	double px = p*std::sin(thetadeg*DEGRAD)*std::cos(phideg*DEGRAD);
	double py = p*std::sin(thetadeg*DEGRAD)*std::sin(phideg*DEGRAD);
	double pz = p*std::cos(thetadeg*DEGRAD);
	return TLorentzVector(px,py,pz,EMeV+massMeV);
}


//helper function to get 3x3 covariance matrix in MeV^2/c^2 for (px, py, pz)
std::array<std::array<double, 3>, 3> SABREPID_N2_M2::GetMomentumCovariance(double E, double theta_deg, double phi_deg, double mass, double sigE, double sigTheta_deg, double sigPhi_deg) const {
	double theta = theta_deg * DEGRAD;
	double phi = phi_deg * DEGRAD;
	double sigma_theta = sigTheta_deg * DEGRAD;
	double sigma_phi = sigPhi_deg * DEGRAD;

	double p = std::sqrt(E*(E+2.*mass));
	if(p <= 0.) return {{{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}}};
	double dP_dE = (E + mass) / p;

	//compute 3x3 jacobian J
	//column 0: d/dE
	//column 1: d/dTheta
	//column 2: d/dPhi
	double J[3][3] = {
		{dP_dE * std::sin(theta) * std::cos(phi), p * std::cos(theta) * std::cos(phi), -p * std::sin(theta) * std::sin(phi)},
		{dP_dE * std::sin(theta) * std::sin(phi), p * std::cos(theta) * std::sin(phi),  p * std::sin(theta) * std::cos(phi)},
		{dP_dE * std::cos(theta),                -p * std::sin(theta),                  0.}
	};

	//V matrix, diagonal elements:
	double var_E = sigE * sigE;
	double var_theta = sigma_theta * sigma_theta;
	double var_phi = sigma_phi * sigma_phi;

	//compute sigma_p = J * V * J^T
	std::array<std::array<double, 3>, 3> Cov = {{{0}}};
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			Cov[i][j] = J[i][0] * var_E * J[j][0] +
						J[i][1] * var_theta * J[j][1] +
						J[i][2] * var_phi * J[j][2];
		}
	}

	return Cov;
}

//correlated 3momentum residual chi^2:
double SABREPID_N2_M2::ComputeCovarianceChi2(const std::array<std::array<double,3>,3>& Cov, const TLorentzVector& residual) const {

	///symmetric covariance matrix:
	//
	// |a b c|
	// |b d e|
	// |c e f|

	double a = Cov[0][0];
	double b = Cov[0][1];
	double c = Cov[0][2];
	double d = Cov[1][1];
	double e = Cov[1][2];
	double f = Cov[2][2];

	//determinant:
	double det = a*(d*f - e*e) - b*(b*f - c*e) + c*(b*e - c*d);
	if(std::abs(det) < 1e-20){
		return 1e9;
	}

	//inverse covariance matrix elements
	double inv00 = (d*f - e*e) / det;
	double inv01 = (c*e - b*f) / det;
	double inv02 = (b*e - c*d) / det;

	double inv11 = (a*f - c*c) / det;
	double inv12 = (b*c - a*e) / det;

	double inv22 = (a*d - b*b) / det;

	double x = residual.Px();
	double y = residual.Py();
	double z = residual.Pz();

	// r^T Sigma^-1 r
	double chi2 = x*(inv00*x + inv01*y + inv02*z)
				+ y*(inv01*x + inv11*y + inv12*z)
				+ z*(inv02*x + inv12*y + inv22*z);

	//return (chi2 > 0.0) ? chi2 : 1e9;
	return (std::isfinite(chi2) && chi2 >= 0.) ? chi2 : 1e9;

}

//evalute single assignment:
double SABREPID_N2_M2::ComputePermutationChi2(const std::array<int, 2>& perm, const double E[2], const double theta[2], const double phi[2], double SPS_E, double SPSTheta, double SPSPhi, const TLorentzVector& P_recoil_expected, TLorentzVector& P_residual_out) const {
	std::array<std::array<double, 3>, 3> Sigma_tot{};
	TLorentzVector P_decay_sum(0.0, 0.0, 0.0, 0.0);

	// Add the two reconstructed decay-product four-vectors and covariance
	// contributions.
	for (int i = 0; i < 2; ++i) {
		const int hitIndex = perm[i];
		const double mass = fHypothesis.final_masses[i];

		const TLorentzVector P_particle = Build4Vector(
			E[hitIndex],
			theta[hitIndex],
			phi[hitIndex],
			mass
		);

		P_decay_sum += P_particle;

		const auto Cov_particle = GetMomentumCovariance(
			E[hitIndex],
			theta[hitIndex],
			phi[hitIndex],
			mass,
			fSigmaE,
			fSigmaTheta,
			fSigmaPhi
		);

		for (int r = 0; r < 3; ++r) {
			for (int c = 0; c < 3; ++c) {
				Sigma_tot[r][c] += Cov_particle[r][c];
			}
		}
	}

	// the expected recoil is reconstructed using the measured SPS ejectile
	// its associated measurement uncertainty therefore contributes to the residual covariance matrix!
	const auto Cov_SPS = GetMomentumCovariance(SPS_E, SPSTheta, SPSPhi, fHypothesis.mass_ejectile, fSPSSigmaE, fSPSSigmaTheta, fSPSSigmaPhi);

	for (int r = 0; r < 3; ++r) {
		for (int c = 0; c < 3; ++c) {
			Sigma_tot[r][c] += Cov_SPS[r][c];
		}
	}

	P_residual_out = P_recoil_expected - P_decay_sum;

	return ComputeCovarianceChi2(Sigma_tot, P_residual_out);
}

void SABREPID_N2_M2::InitDiagnostics(TDirectory* targetDir) {
	if (!targetDir) {
		return;
	}

	outdir = targetDir->mkdir("SABREPID_N2_M2_Diagnostics");

	if (!outdir) outdir = targetDir->GetDirectory("SABREPID_N2_M2_Diagnostics");

	if (!outdir) return;

	outdir->cd();

	hBestChi2 = new TH1D("hBestChi2","Best assignment reduced #chi^{2};#chi^{2}/3;Counts",500,0.0,20.0);


	hBestPermutation = new TH1D("hBestPermutation","Best particle assignment;Assignment;Counts",2,-0.5,1.5);
	hChi2_BestVsNext = new TH2D("hChi2_BestVsNext","Best #chi^{2} vs next-best #chi^{2};Best #chi^{2}/3;Next-best #chi^{2}/3",500,0.0,20.0,500,0.0,20.0);

	h2Chi2ByPermutation = new TH2D("h2Chi2ByPermutation", "#chi^{2} by particle assignment;Assignment;#chi^{2}", 2, -0.5, 1.5, 500, 0.0, 500.0);

	h2Chi2DifByPermutation = new TH2D("h2Chi2DifByPermutation","(#chi^{2}_{i}-#chi^{2}_{best}) by assignment;Assignment;#Delta#chi^{2}",2,-0.5,1.5,500,0.0,500.0);

	for(int i = 0; i < 2; ++i) {
		const TString label = Form("P%d: (%s,%s)", i, fHypothesis.final_particles[allPerms[i][0]].Data(), fHypothesis.final_particles[allPerms[i][1]].Data());

		hBestPermutation->GetXaxis()->SetBinLabel(i + 1, label.Data());
		h2Chi2ByPermutation->GetXaxis()->SetBinLabel(i + 1, label.Data());
		h2Chi2DifByPermutation->GetXaxis()->SetBinLabel(i + 1, label.Data());
	}
}

PIDResult_N2_M2 SABREPID_N2_M2::EvaluateEvent(const double E[2],const double theta[2],const double phi[2],double SPS_E,double SPSTheta,double SPSPhi) {
	PIDResult_N2_M2 res;
	res.Reset();

	// assuming beam aligned with +z.
	double beam_p = std::sqrt(fHypothesis.beamEnergyMeV*(fHypothesis.beamEnergyMeV + 2.*fHypothesis.mass_beam));

	TLorentzVector P_beam(0.,0.,beam_p,fHypothesis.beamEnergyMeV + fHypothesis.mass_beam);

	TLorentzVector P_target(0.,0.,0.,fHypothesis.mass_target);

	TLorentzVector P_ejectile = Build4Vector(SPS_E,SPSTheta,SPSPhi,fHypothesis.mass_ejectile);

	//reconstructed four momentum of the recoil populated in reactin before its 2 body decay
	TLorentzVector P_recoil_expected = P_beam + P_target - P_ejectile;

	double min_chi2 = 1e9;
	double second_best_chi2 = 1e9;

	int best_perm_index = -1;
	TLorentzVector best_residual;

	// P0 = particle 0 assigned to SABRE hit 0, particle 1 to hit 1
	// P1 = particle 0 assigned to SABRE hit 1, particle 1 to hit 0
	for (size_t permIndex = 0; permIndex < allPerms.size(); permIndex++) {
		const auto& perm = allPerms[permIndex];

		TLorentzVector P_residual;

		double chi2 = ComputePermutationChi2(perm,E,theta,phi,SPS_E,SPSTheta,SPSPhi,P_recoil_expected,P_residual);

		res.permChi2s[permIndex] = chi2;

		if (h2Chi2ByPermutation) h2Chi2ByPermutation->Fill(static_cast<double>(permIndex),chi2);

		if (chi2 < min_chi2) {
			second_best_chi2 = min_chi2;
			min_chi2 = chi2;
			best_perm_index = static_cast<int>(permIndex);
			best_residual = P_residual;
		} else if(chi2 < second_best_chi2) {
			second_best_chi2 = chi2;
		}
	}

	if(best_perm_index < 0) return res;

	res.bestChi2Index = best_perm_index;
	res.bestChi2 = min_chi2;
	res.hit_indices = allPerms[best_perm_index];

	res.missing_px = best_residual.Px();
	res.missing_pz = best_residual.Pz();
	res.missing_py = best_residual.Py();
	res.missing_E = best_residual.E();
	res.missing_Pmag = best_residual.P();

	const double reducedChi2 = res.bestChi2 / 3.0;
	res.passesCut = (reducedChi2 <= fMaxChi2Cut);

	if(hBestPermutation) hBestPermutation->Fill(best_perm_index);

	if(hBestChi2) hBestChi2->Fill(reducedChi2);

	if(hChi2_BestVsNext && second_best_chi2 < 1e9) hChi2_BestVsNext->Fill(reducedChi2,second_best_chi2 / 3.0);

	//fill only after determing best chi2 (so after cycling all perms once)
	if (h2Chi2DifByPermutation) {
		for(size_t permIndex = 0; permIndex < allPerms.size(); permIndex++) {
			double deltaChi2 = res.permChi2s[permIndex] - res.bestChi2;
			h2Chi2DifByPermutation->Fill(static_cast<double>(permIndex),deltaChi2);
		}
	}

	return res;
}
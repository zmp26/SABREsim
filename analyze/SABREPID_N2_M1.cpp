#include "SABREPID_N2_M1.h"
#include "PID_structs.h"

#include <cmath>

#include "TString.h"

SABREPID_N2_M1::SABREPID_N2_M1() = default;

SABREPID_N2_M1::~SABREPID_N2_M1() = default;

void SABREPID_N2_M1::SetHypothesis(const PIDHypothesis_N2& hypo){
	fHypothesis = hypo;
}

void SABREPID_N2_M1::SetResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg){
	fSigmaE = sigEMeV;
	fSigmaTheta = sigThetaDeg;
	fSigmaPhi = sigPhiDeg;
}

void SABREPID_N2_M1::SetSPSResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg){
	fSPSSigmaE = sigEMeV;
	fSPSSigmaTheta = sigThetaDeg;
	fSPSSigmaPhi = sigPhiDeg;
}

TLorentzVector SABREPID_N2_M1::Build4Vector(double EMeV, double thetaDeg, double phideg, double massMeV) const {
	double p = std::sqrt(EMeV*EMeV + 2.*EMeV*massMeV);
	double px = p*std::sin(thetaDeg*DEGRAD)*std::cos(phideg*DEGRAD);
	double py = p*std::sin(thetaDeg*DEGRAD)*std::sin(phideg*DEGRAD);
	double pz = p*std::cos(thetaDeg*DEGRAD);
	return TLorentzVector(px, py, pz, EMeV+massMeV);
}

//computes a 4x4 covariance matrix in MeV^2 for (px, py, pz, E_tot)
std::array<std::array<double, 4>, 4> SABREPID_N2_M1::Get4MomentumCovariance(
	double E, double theta_deg, double phi_deg, double mass, 
	double sigE, double sigTheta_deg, double sigPhi_deg) const {

	double theta = theta_deg * DEGRAD;
	double phi = phi_deg * DEGRAD;
	double sigma_theta = sigTheta_deg * DEGRAD;
	double sigma_phi = sigPhi_deg * DEGRAD;

	double p = std::sqrt(E * (E + 2. * mass));
	if (p <= 0.) return {{{0}}};

	double dP_dE = (E + mass) / p;

	//4x3 Jacobian J = d(px, py, pz, E_tot) / d(E_kin, theta, phi)
	double J[4][3] = {
		{dP_dE * std::sin(theta) * std::cos(phi),	p * std::cos(theta) * std::cos(phi), -p * std::sin(theta) * std::sin(phi)},
		{dP_dE * std::sin(theta) * std::sin(phi),	p * std::cos(theta) * std::sin(phi),  p * std::sin(theta) * std::cos(phi)},
		{dP_dE * std::cos(theta), -p * std::sin(theta), 0.0},
		{1.0, 0.0, 0.0}
	};

	double var_E = sigE * sigE;
	double var_theta = sigma_theta * sigma_theta;
	double var_phi = sigma_phi * sigma_phi;

	std::array<std::array<double, 4>, 4> Cov = {{{0}}};
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			Cov[i][j] = J[i][0] * var_E * J[j][0] + J[i][1] * var_theta * J[j][1] + J[i][2] * var_phi * J[j][2];
		}
	}

	return Cov;
}

double SABREPID_N2_M1::ComputePermutationChi2(
	const std::array<int, 2>& perm, const double E[1], const double theta[1], const double phi[1], 
	double SPS_E, double SPSTheta, double SPSPhi, 
	TLorentzVector& P_missing_out, double& missingMass_out) const {

	const int detected_spec = perm[0];
	const int missing_spec = perm[1];

	const double detected_mass = fHypothesis.final_masses[detected_spec];
	const double expected_missing_mass = fHypothesis.final_masses[missing_spec];

	double beam_p = std::sqrt(fHypothesis.beamEnergyMeV * (fHypothesis.beamEnergyMeV + 2. * fHypothesis.mass_beam));
	TLorentzVector P_beam(0., 0., beam_p, fHypothesis.beamEnergyMeV + fHypothesis.mass_beam);
	TLorentzVector P_target(0., 0., 0., fHypothesis.mass_target);

	TLorentzVector P_ejectile = Build4Vector(SPS_E, SPSTheta, SPSPhi, fHypothesis.mass_ejectile);
	TLorentzVector P_detected = Build4Vector(E[0], theta[0], phi[0], detected_mass);

	// Reconstruct missing particle 4-momentum
	P_missing_out = P_beam + P_target - P_ejectile - P_detected;

	double M2 = P_missing_out.M2();
	if (M2 <= 0.0) {
		missingMass_out = 0.0;
		return 1e9;
	}

	missingMass_out = std::sqrt(M2);

	//sum covariance matrix contributions from SPS ejectile and detected SABRE particle
	auto Cov_SPS = Get4MomentumCovariance(SPS_E, SPSTheta, SPSPhi, fHypothesis.mass_ejectile, fSPSSigmaE, fSPSSigmaTheta, fSPSSigmaPhi);
	auto Cov_det = Get4MomentumCovariance(E[0], theta[0], phi[0], detected_mass, fSigmaE, fSigmaTheta, fSigmaPhi);

	std::array<std::array<double, 4>, 4> Cov_missing{{{0}}};
	for (int r = 0; r < 4; ++r) {
		for (int c = 0; c < 4; ++c) {
			Cov_missing[r][c] = Cov_SPS[r][c] + Cov_det[r][c];
		}
	}

	//gradient vector of missing invariant mass M = sqrt(E_tot^2 - px^2 - py^2 - pz^2)
	double px = P_missing_out.Px();
	double py = P_missing_out.Py();
	double pz = P_missing_out.Pz();
	double Etot = P_missing_out.E();

	std::array<double, 4> grad = { -px / missingMass_out, -py / missingMass_out, -pz / missingMass_out, Etot / missingMass_out };

	//transform 4momentum covariance to scalar mass variance: sigma_M^2 = grad^T * Cov_missing * grad
	double var_M = 0.0;
	for (int r = 0; r < 4; ++r) {
		for (int c = 0; c < 4; ++c) {
			var_M += grad[r] * Cov_missing[r][c] * grad[c];
		}
	}

	if (!std::isfinite(var_M) || var_M <= 0.0) return 1e9;

	double diff = missingMass_out - expected_missing_mass;
	double chi2 = (diff * diff) / var_M;

	return (std::isfinite(chi2) && chi2 >= 0.0) ? chi2 : 1e9;
}

void SABREPID_N2_M1::InitDiagnostics(TDirectory* targetDir) {
	if (!targetDir) return;

	outdir = targetDir->mkdir("SABREPID_N2_M1_Diagnostics");
	if (!outdir) outdir = targetDir->GetDirectory("SABREPID_N2_M1_Diagnostics");
	if (!outdir) return;

	outdir->cd();

	hBestChi2 = new TH1D("hBestChi2", "Best assignment reduced #chi^{2};#chi^{2};Counts", 500, 0.0, 20.0);
	hBestPermutation = new TH1D("hBestPermutation", "Best particle assignment;Assignment;Counts", 2, -0.5, 1.5);
	hChi2_BestVsNext = new TH2D("hChi2_BestVsNext", "Best #chi^{2} vs next-best #chi^{2};Best #chi^{2};Next-best #chi^{2}", 500, 0.0, 20.0, 500, 0.0, 20.0);
	h2Chi2ByPermutation = new TH2D("h2Chi2ByPermutation", "#chi^{2} by particle assignment;Assignment;#chi^{2}", 2, -0.5, 1.5, 500, 0.0, 500.0);
	h2Chi2DifByPermutation = new TH2D("h2Chi2DifByPermutation", "(#chi^{2}_{i}-#chi^{2}_{best}) by assignment;Assignment;#Delta#chi^{2}", 2, -0.5, 1.5, 500, 0.0, 500.0);
	hMissingMass = new TH1D("hMissingMass", "Reconstructed Missing Mass;Mass (MeV/c^{2});Counts", 1000, 0.0, 10000.0);

	for (size_t i = 0; i < allPerms.size(); ++i) {
		const TString label = Form("Det: %s, Miss: %s", fHypothesis.final_particles[allPerms[i][0]].Data(), fHypothesis.final_particles[allPerms[i][1]].Data());
		hBestPermutation->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, label.Data());
		h2Chi2ByPermutation->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, label.Data());
		h2Chi2DifByPermutation->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, label.Data());
	}
}

PIDResult_N2_M1 SABREPID_N2_M1::EvaluateEvent(const double E[1], const double theta[1], const double phi[1], double SPS_E, double SPSTheta, double SPSPhi) {
	PIDResult_N2_M1 res;
	res.Reset();

	double min_chi2 = 1e9;
	double second_best_chi2 = 1e9;
	int best_perm_index = -1;

	TLorentzVector best_P_missing;
	double best_missing_mass = 0.0;

	for (size_t permIndex = 0; permIndex < allPerms.size(); ++permIndex) {
		const auto& perm = allPerms[permIndex];
		TLorentzVector P_missing;
		double missing_mass = 0.0;

		double chi2 = ComputePermutationChi2(perm, E, theta, phi, SPS_E, SPSTheta, SPSPhi, P_missing, missing_mass);

		res.permChi2s[permIndex] = chi2;

		if (h2Chi2ByPermutation) h2Chi2ByPermutation->Fill(static_cast<double>(permIndex), chi2);

		if (chi2 < min_chi2) {
			second_best_chi2 = min_chi2;
			min_chi2 = chi2;
			best_perm_index = static_cast<int>(permIndex);
			best_P_missing = P_missing;
			best_missing_mass = missing_mass;
		} else if (chi2 < second_best_chi2) {
			second_best_chi2 = chi2;
		}
	}

	if (best_perm_index < 0) return res;

	res.bestChi2Index = best_perm_index;
	res.bestChi2 = min_chi2;
	res.detected_hit_index = 0;
	res.detected_species_index = allPerms[best_perm_index][0];
	res.missing_species_index = allPerms[best_perm_index][1];

	res.missing_px = best_P_missing.Px();
	res.missing_py = best_P_missing.Py();
	res.missing_pz = best_P_missing.Pz();
	res.missing_E = best_P_missing.E();
	res.missing_Pmag = best_P_missing.P();
	res.missing_MassCalc = best_missing_mass;

	//1 degree of freedom (invariant mass constraint)
	const double reducedChi2 = res.bestChi2;
	res.passesCut = (reducedChi2 <= fMaxChi2Cut);

	if(hBestPermutation) hBestPermutation->Fill(best_perm_index);
	if(hBestChi2) hBestChi2->Fill(reducedChi2);
	if(hMissingMass) hMissingMass->Fill(best_missing_mass);
	if(hChi2_BestVsNext && second_best_chi2 < 1e9) hChi2_BestVsNext->Fill(reducedChi2, second_best_chi2);

	if(h2Chi2DifByPermutation) {
		for(size_t permIndex = 0; permIndex < allPerms.size(); ++permIndex) {
			double deltaChi2 = res.permChi2s[permIndex] - res.bestChi2;
			h2Chi2DifByPermutation->Fill(static_cast<double>(permIndex), deltaChi2);
		}
	}

	return res;
}
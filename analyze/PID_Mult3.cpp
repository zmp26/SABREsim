#include "PID_Mult3.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>

constexpr double DEGRAD = M_PI/180.;
constexpr double RADDEG = 180./M_PI;

//helper function to get 3x3 covariance matrix in MeV^2/c^2 for (px, py, pz)
std::array<std::array<double, 3>, 3> PID_Mult3::GetMomentumCovariance(double E, double theta_deg, double phi_deg, double mass, double sigE, double sigTheta_deg, double sigPhi_deg) const {
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

double PID_Mult3::ComputeCovarianceChi2(const std::array<std::array<double,3>,3>& Cov, const TLorentzVector& residual) const {

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

	return (chi2 > 0.0) ? chi2 : 1e9;

}

PID_Mult3::PID_Mult3(){
	//constructor
}

PID_Mult3::~PID_Mult3(){
	//destructor
	//delete hBestChi2;
	//delete hChi2_BestVsNext;
}

TLorentzVector PID_Mult3::Build4Vector(double EMeV, double thetadeg, double phideg, double massMeV) const {
	double p = std::sqrt(EMeV*EMeV + 2.*EMeV*massMeV);//p = sqrt(E^2 + 2*E*M)
	double px = p*std::sin(thetadeg*DEGRAD)*std::cos(phideg*DEGRAD);
	double py = p*std::sin(thetadeg*DEGRAD)*std::sin(phideg*DEGRAD);
	double pz = p*std::cos(thetadeg*DEGRAD);
	return TLorentzVector(px, py, pz, EMeV+massMeV);
}

void PID_Mult3::SetHypothesis(const PIDHypothesis& hypo){
	fHypothesis = hypo;
}

void PID_Mult3::SetResolution(double sigE, double sigThetaDeg, double sigPhiDeg){
	fSigmaE = sigE;
	fSigmaTheta = sigThetaDeg;
	fSigmaPhi = sigPhiDeg;
}

void PID_Mult3::SetSPSResolution(double sigE, double sigTheta, double sigPhi){
	fSPSSigmaE = sigE;
	fSPSSigmaTheta = sigTheta;
	fSPSSigmaPhi = sigPhi;
}

void PID_Mult3::InitDiagnostics(TDirectory* targetDir){
	if(!targetDir) return;

	outdir = targetDir->mkdir("PID_Mult3_Diagnostics");
	if(!outdir){
		outdir = targetDir->GetDirectory("PID_Mult3_Diagnostics");
	}

	outdir->cd();
	hBestChi2 = new TH1D("hBestChi2", "Best Permutation Reduced #chi^{2}", 500, 0, 5);
	hBestPermutation = new TH1D("hBestPermutation", "Permutation with lowest #chi^{2};Permutation;Counts", 6, -0.5, 5.5);
	hChi2_BestVsNext = new TH2D("hChi2_BestVsNext", "Best #chi^{2} vs 2nd Best #chi^{2};Best #chi^{2};2nd Best #chi^{2}", 500, 0, 1, 500, 0, 1);
	h2Chi2ByPermutation = new TH2D("h2Chi2ByPermutation", "#chi^{2} by permutation;permutation;#chi^{2}", 6, -0.5, 5.5, 500, 0, 500);

	for(int i=0; i<6; i++){
		TString label = Form("P%d: (%s,%s,%s)", i, fHypothesis.final_particles[allPerms[i][0]].Data(), fHypothesis.final_particles[allPerms[i][1]].Data(), fHypothesis.final_particles[allPerms[i][2]].Data());
		h2Chi2ByPermutation->GetXaxis()->SetBinLabel(i+1, label.Data());
		hBestPermutation->GetXaxis()->SetBinLabel(i+1, label.Data());
	}
}

double PID_Mult3::ComputePermutationChi2(const std::array<int,3>& perm,
								  const double E[3], const double theta[3], const double phi[3],
								  const double SPS_E, const double SPSTheta, const double SPSPhi,
								  const TLorentzVector& P_recoil_expected,
								  TLorentzVector& P_residual_out) const {

	//perm[0] = hit index for final_masses[0]
	//perm[1] = hit index for final_masses[1]
	//perm[2] = hit index for final_masses[2]
	std::array<std::array<double,3>,3> Sigma_tot{};
	TLorentzVector P_decay_sum(0., 0., 0., 0.);
	double var_px = 0.;
	double var_py = 0.;
	double var_pz = 0.;
	double var_E = 0.;

	for(int i=0; i<3; i++){
		int hitindex = perm[i];
		double mass = fHypothesis.final_masses[i];

		TLorentzVector p4 = Build4Vector(E[hitindex], theta[hitindex], phi[hitindex], mass);
		P_decay_sum += p4;

		//fetch pixel-dependent resolutions here (eventually, we add pixel handler here so we have a map ready!)
		double sigE = fSigmaE;//		 use constant class value for now
		double sigTheta = fSigmaTheta;//"								 "
		double sigPhi = fSigmaPhi;//    "								 "

		//add particle covariance matrix to total covariance matrix
		auto Cov_particle = GetMomentumCovariance(E[hitindex], theta[hitindex], phi[hitindex], mass, sigE, sigTheta, sigPhi);
		for(int r=0; r<3; r++){
			for(int c=0; c<3; c++){
				Sigma_tot[r][c] += Cov_particle[r][c];
			}
		}
	}

	//add SPS ejectile covariance matrix
	double sps_sigE = fSPSSigmaE;
	double sps_sigTheta = fSPSSigmaTheta;
	double sps_sigPhi = fSPSSigmaPhi;

	auto Cov_SPS = GetMomentumCovariance(SPS_E, SPSTheta, SPSPhi, fHypothesis.mass_ejectile, sps_sigE, sps_sigTheta, sps_sigPhi);
	for(int r=0; r<3; r++){
		for(int c=0; c<3; c++){
			Sigma_tot[r][c] += Cov_SPS[r][c];
		}
	}

	//compute missing momentum
	P_residual_out = P_recoil_expected - P_decay_sum;

/*
	//grab diagonals for individual Px, Py, Pz variances:
	var_px = Sigma_tot[0][0];
	var_py = Sigma_tot[1][1];
	var_pz = Sigma_tot[2][2];

	double chi2_px = (P_residual_out.Px() * P_residual_out.Px()) / (var_px > 0. ? var_px : 1.);
	double chi2_py = (P_residual_out.Py() * P_residual_out.Py()) / (var_py > 0. ? var_py : 1.);
	double chi2_pz = (P_residual_out.Pz() * P_residual_out.Pz()) / (var_pz > 0. ? var_pz : 1.);

	//Energy residula contribution
	//var_E = 3.*fSigmaE*fSigmaE + fSPSSigmaE*fSPSSigmaE;
	//double chi2_E = (P_residual_out.E() * P_residual_out.E()) / (var_E > 0. ? var_E : 1.);

	return (chi2_px + chi2_py + chi2_pz);// + chi2_E);
*/

	return ComputeCovarianceChi2(Sigma_tot, P_residual_out);
}

PIDResult PID_Mult3::EvaluateEvent(
	const double E[3], const double theta[3], const double phi[3],
	double SPS_E, double SPSTheta, double SPSPhi) 
{
	PIDResult res;
	res.Reset();

	// Reconstruct expected 4-momentum of recoil from initial kinematics
	double beam_p = std::sqrt(fHypothesis.beamEnergyMeV * (fHypothesis.beamEnergyMeV + 2.0 * fHypothesis.mass_beam));
	TLorentzVector P_beam(0.0, 0.0, beam_p, fHypothesis.beamEnergyMeV + fHypothesis.mass_beam);
	TLorentzVector P_target(0.0, 0.0, 0.0, fHypothesis.mass_target);
	TLorentzVector P_ejectile = Build4Vector(SPS_E, SPSTheta, SPSPhi, fHypothesis.mass_ejectile);

	TLorentzVector P_recoil_expected = P_beam + P_target - P_ejectile;

	double min_chi2 = 1e9;
	double second_best_chi2 = 1e9;
	int best_perm_index = -1;
	TLorentzVector best_p_missing;

	// Evaluate all 6 permutations
	for (size_t permindex = 0; permindex < allPerms.size(); ++permindex) {
		const auto& perm = allPerms[permindex];

		TLorentzVector P_missing;
		double chi2 = ComputePermutationChi2(
			perm, E, theta, phi, 
			SPS_E, SPSTheta, SPSPhi, 
			P_recoil_expected, P_missing
		);

		// Store Chi2 per permutation
		res.permChi2s[permindex] = chi2;
		h2Chi2ByPermutation->Fill(permindex, chi2);

		if (chi2 < min_chi2) {
			second_best_chi2 = min_chi2;
			min_chi2 = chi2;
			best_perm_index = static_cast<int>(permindex);
			best_p_missing = P_missing;
		} else if (chi2 < second_best_chi2) {
			second_best_chi2 = chi2;
		}
	}

	if (best_perm_index >= 0) {
		hBestPermutation->Fill(best_perm_index);
		res.bestChi2Index = best_perm_index;
		res.bestChi2 = min_chi2;
		res.hit_indices = allPerms[best_perm_index];

		res.missing_px   = best_p_missing.Px();
		res.missing_py   = best_p_missing.Py();
		res.missing_pz   = best_p_missing.Pz();
		res.missing_E    = best_p_missing.E();
		res.missing_Pmag = best_p_missing.P();

		// Reduced Chi2 with 4 degrees of freedom (px, py, pz)
		double reducedChi2 = res.bestChi2 / 3.0;
		res.passesCut = (reducedChi2 <= fMaxChi2Cut);

		if (outdir && hBestChi2) {
			hBestChi2->Fill(reducedChi2);
			if (hChi2_BestVsNext && second_best_chi2 < 1e8) {
				hChi2_BestVsNext->Fill(reducedChi2, second_best_chi2 / 3.0);
			}
		}
	}

	return res;
}
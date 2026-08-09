#include "PID_Mult2.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>

PID_Mult2::PID_Mult2(){
	//constructor
}

PID_Mult2::~PID_Mult2(){
	//destructor
}

TLorentzVector PID_Mult2::Build4Vector(double EMeV, double thetadeg, double phideg, double massMeV) const {
	double p = std::sqrt(EMeV*EMeV + 2.*EMeV*massMeV);
	double px = p*std::sin(thetadeg*DEGRAD)*std::cos(phideg*DEGRAD);
	double py = p*std::sin(thetadeg*DEGRAD)*std::sin(phideg*DEGRAD);
	double pz = p*std::cos(thetadeg*DEGRAD);
	return TLorentzVector(px,py,pz,EMeV+massMeV);
}

std::array<std::array<double,4>,4> PID_Mult2::GetP4Covariance(double EMeV, double thetadeg, double phideg, double massMeV,
															  double sigE, double sigThetaDeg, double sigPhiDeg) const {

	double theta = thetadeg * DEGRAD;
	double phi = phideg*DEGRAD;
	double sigtheta = sigThetaDeg*DEGRAD;
	double sigphi = sigPhiDeg*DEGRAD;

	double p = std::sqrt(EMeV*EMeV + 2.*EMeV*massMeV);
	if(p<=0.) return {{{0.}}};

	double dP_dE = (EMeV+massMeV)/p;//dP/dE = (E + M)/sqrt(E*E + 2.*E*M) = (E+M)/p

	//jacobian matrix is 4x3: rows = (E, px, py, pz). cols = (E_meas, theta_meas, phi_meas)
	double J[4][3] = {
		{1.0, 0., 0.},//dE/dE = 1, dE/dTheta = 0, dE/dPhi = 0
		{dP_dE*std::sin(theta)*std::cos(phi), p*std::cos(theta)*std::cos(phi), -p*std::sin(theta)*std::sin(phi)},
		{dP_dE*std::sin(theta)*std::sin(phi), p*std::cos(theta)*std::sin(phi), p*std::sin(theta)*std::cos(phi)},
		{dP_dE*std::cos(theta),              -p*std::sin(theta),               0.}
	};

	double var[3] = {sigE*sigE, sigtheta*sigtheta, sigphi*sigphi};

	//Cov = J * V * J^T
	std::array<std::array<double,4>,4> Cov = {{{0}}};
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			for(int k=0; k<3; k++){
				Cov[i][j] += J[i][k] * var[k] * J[j][k];
			}
		}
	}
	return Cov;
}

void PID_Mult2::SetHypothesis(const PIDHypothesis& hypo){
	fHypothesis = hypo;
}

void PID_Mult2::SetResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg){
	fSigmaE = sigEMeV;
	fSigmaTheta = sigThetaDeg;
	fSigmaPhi = sigPhiDeg;
}

void PID_Mult2::SetSPSResolution(double sigEMeV, double sigThetaDeg, double sigPhiDeg){
	fSPSSigmaE = sigEMeV;
	fSPSSigmaTheta = sigThetaDeg;
	fSPSSigmaPhi = sigPhiDeg;
}

void PID_Mult2::InitDiagnostics(TDirectory* targetDir){
	if(!targetDir) return;

	outdir = targetDir->mkdir("PID_Mult2_Diagnostics");
	if(!outdir) outdir = targetDir->GetDirectory("PID_Mult2_Diagnostics");

	outdir->cd();
	hBestChi2 = new TH1D("hBestChi2", "Best Permutation #chi^{2};#chi^{2};Counts", 500, 0, 50);
	hBestPermutation = new TH1D("hBestPermutation", "Permutation with lowest #chi^{2};Permutation;Counts", 6, -0.5, 5.5);
	hMissingMassBest = new TH1D("hMissingMassBest", "Reconstructed Missing Mass (Best Permutation);M_{calc} [MeV/c^{2}];Counts", 500, 0, 3000);
	hSigmaMBest = new TH1D("hSigmaMBest", "Calculated Mass Resolution #sigma_{M} (Best Permutation);#sigma_{M} [MeV/c^{2}];Counts", 200, 0, 20);
	hChi2_BestVsNext = new TH2D("hChi2_BestVsNext", "Best #chi^{2} vs 2nd Best #chi^{2};Best #chi^{2};2nd Best #chi^{2}", 500, 0, 50, 500, 0, 50);
	h2Chi2ByPermutation = new TH2D("h2Chi2ByPermutation", "#chi^{2} by permutation;Permutation;#chi^{2}", 6, -0.5, 5.5, 500, 0, 500);

	for(int i=0; i<6; i++){
		TString label = Form("P%d: H0=%s, H1=%s, Miss=%s", i, fHypothesis.final_particles[allPerms[i][0]].Data(), fHypothesis.final_particles[allPerms[i][1]].Data(),fHypothesis.final_particles[allPerms[i][2]].Data());
		h2Chi2ByPermutation->GetXaxis()->SetBinLabel(i+1, label.Data());
		hBestPermutation->GetXaxis()->SetBinLabel(i+1, label.Data());
	}
}

double PID_Mult2::ComputePermutationChi2(const std::array<int,3>& perm, const double E[2], const double theta[2], const double phi[2],
										 double SPS_E, double SPSTheta, double SPSPhi,
										 const TLorentzVector& P_recoil_expected,
										 TLorentzVector& P_missing_out,
										 double& sigma_M_out) const {

	//perm[0] = species index for Detected Hit 0 (H0)
	//perm[1] = species index for Detected Hit 1 (H1)
	//perm[2] = species index for Missing Particle (M)

	double mass_H0 = fHypothesis.final_masses[perm[0]];
	double mass_H1 = fHypothesis.final_masses[perm[1]];
	double mass_M_nominal = fHypothesis.final_masses[perm[2]];

	TLorentzVector P4_H0 = Build4Vector(E[0], theta[0], phi[0], mass_H0);
	TLorentzVector P4_H1 = Build4Vector(E[1], theta[1], phi[1], mass_H1);

	P_missing_out = P_recoil_expected - P4_H0 - P4_H1;

	double M2_calc = P_missing_out.M2();
	double M_calc = P_missing_out.M();
	if(M_calc < 0. || std::isnan(M_calc)){
		sigma_M_out = 1.;
		return 1e9;
	}

	//calculate the combined covariance matrix:
	auto Cov_H0 = GetP4Covariance(E[0], theta[0], phi[0], mass_H0, fSigmaE, fSigmaTheta, fSigmaPhi);
	auto Cov_H1 = GetP4Covariance(E[1], theta[1], phi[1], mass_H1, fSigmaE, fSigmaTheta, fSigmaPhi);
	auto Cov_SPS = GetP4Covariance(SPS_E, SPSTheta, SPSPhi, fHypothesis.mass_ejectile, fSPSSigmaE, fSPSSigmaTheta, fSPSSigmaPhi);

	std::array<std::array<double,4>,4> Sigma_P4_tot = {{{0}}};
	for(int r=0; r<4; r++){
		for(int c=0; c<4; c++){
			Sigma_P4_tot[r][c] = Cov_H0[r][c] + Cov_H1[r][c] + Cov_SPS[r][c];
		}
	}

	//compute gradient vector g = DM/dP4 = (E_miss, -px_miss, -py_miss, -pz_miss) / M_calc
	double g[4] = {
		P_missing_out.E()  / M_calc,
		-P_missing_out.Px() / M_calc,
		-P_missing_out.Py() / M_calc,
		-P_missing_out.Pz() / M_calc
	};

	//propagate error: var_M = g^T * Sigma_P4_tot * g
	double var_M = 0.;
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			var_M += g[i] * Sigma_P4_tot[i][j] * g[j];
		}
	}


	sigma_M_out = (var_M > 0.) ? std::sqrt(var_M) : 1.0;

	//evaluate mass Chi2 with event-specific resolution
	double mass_diff = M_calc - mass_M_nominal;
	double chi2_mass = (mass_diff * mass_diff) / (sigma_M_out * sigma_M_out);

	return chi2_mass;
}

PIDResult_Mult2 PID_Mult2::EvaluateEvent(const double E[2], const double theta[2], const double phi[2],
								  const double SPS_E, const double SPSTheta, const double SPSPhi){

	PIDResult_Mult2 res;
	res.Reset();

	double beam_p = std::sqrt(fHypothesis.beamEnergyMeV * fHypothesis.beamEnergyMeV + 2.*fHypothesis.beamEnergyMeV*fHypothesis.mass_beam);
	TLorentzVector P_beam(0., 0., beam_p, fHypothesis.beamEnergyMeV+fHypothesis.mass_beam);
	TLorentzVector P_target(0., 0., 0., fHypothesis.mass_target);
	TLorentzVector P_ejectile = Build4Vector(SPS_E, SPSTheta, SPSPhi, fHypothesis.mass_ejectile);

	TLorentzVector P_recoil_expected = P_beam + P_target - P_ejectile;

	double min_chi2 = 1e9;
	double second_best_chi2 = 1e9;
	int best_perm_index = -1;
	TLorentzVector best_p_missing;
	double best_sigma_M = 0.;

	for(size_t permindex=0; permindex<6; permindex++){
		const auto& perm = allPerms[permindex];

		TLorentzVector P_missing;
		double sigma_M = 0.;
		double chi2 = ComputePermutationChi2(perm, E, theta, phi, SPS_E, SPSTheta, SPSPhi, P_recoil_expected, P_missing, sigma_M);

		res.permChi2s[permindex] = chi2;
		if(h2Chi2ByPermutation){
			h2Chi2ByPermutation->Fill(permindex, chi2);
		}

		if(chi2 < min_chi2){
			second_best_chi2 = min_chi2;
			min_chi2 = chi2;
			best_perm_index = static_cast<int>(permindex);
			best_p_missing = P_missing;
			best_sigma_M = sigma_M;
		} else if(chi2 < second_best_chi2){
			second_best_chi2 = chi2;
		}
	}

	if(best_perm_index >= 0){
		res.bestChi2Index = best_perm_index;
		res.bestChi2 = min_chi2;

		res.hit_indices = {allPerms[best_perm_index][0], allPerms[best_perm_index][1]};
		res.missing_species_index = allPerms[best_perm_index][2];

		res.missing_px = best_p_missing.Px();
		res.missing_py = best_p_missing.Py();
		res.missing_pz = best_p_missing.Pz();
		res.missing_Etot = best_p_missing.E();
		res.missing_Pmag = best_p_missing.P();
		res.missing_MassCalc = best_p_missing.M();

		res.passesCut = (res.bestChi2 <= fMaxChi2Cut);

		if(outdir){
			if(hBestPermutation) hBestPermutation->Fill(best_perm_index);
			if(hBestChi2) hBestChi2->Fill(res.bestChi2);
			if(hMissingMassBest) hMissingMassBest->Fill(res.missing_MassCalc);
			if(hSigmaMBest) hSigmaMBest->Fill(best_sigma_M);
			if(hChi2_BestVsNext && second_best_chi2 < 1e8) hChi2_BestVsNext->Fill(res.bestChi2, second_best_chi2); 
		}
	}

	return res;
}

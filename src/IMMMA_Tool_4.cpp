#include "IMMMA_Tool_4.h"
#include "structs.h"
#include <algorithm>

IMMMA_Tool_4::IMMMA_Tool_4() = default;

//setters
void IMMMA_Tool_4::SetBeamNucleus(int A, TString sym, double mass){
	beam = {A, sym, mass};
}

void IMMMA_Tool_4::SetTargetNucleus(int A, TString sym, double mass){
	target = {A, sym, mass};
}

void IMMMA_Tool_4::SetEjectileNucleus(int A, TString sym, double mass){
	ejectile = {A, sym, mass};
}

void IMMMA_Tool_4::SetRecoilNucleus(int A, TString sym, double mass){
	recoil = {A, sym, mass};
}

void IMMMA_Tool_4::SetBreakupNuclei(const std::vector<Nucleus>& b){
	breakups = b;
}

void IMMMA_Tool_4::SetBeamEnergyMeV(double energy){
	beam_energy = energy;
}

void IMMMA_Tool_4::SetRecoilExEMeV(double energy){
	recoil_ExE = energy;
}

//lab frame vector builders
TLorentzVector IMMMA_Tool_4::BuildBeamLV() const {
	IMMMA_Fragment f{beam.massMeV, beam_energy, 0, 0, false};
	return BuildLab4Vector(f);
}

TLorentzVector IMMMA_Tool_4::BuildTargetLV() const {
	IMMMA_Fragment f{target.massMeV, 0, 0, 0, false};
	return BuildLab4Vector(f);
}

TLorentzVector IMMMA_Tool_4::BuildEjectileLV(double E, double th, double ph) const {
	IMMMA_Fragment f{ejectile.massMeV, E, th, ph, false};
	return BuildLab4Vector(f);
}

IMMMA_DecayResult4 IMMMA_Tool_4::SolveSequentialDecay(const TLorentzVector& recoil,
													  const IMMMA_Fragment& fbu1,
													  const IMMMA_Fragment& fbu2,
													  const IMMMA_Fragment& fbu3,
													  TLorentzVector& outLV_bu1,
													  TLorentzVector& outLV_bu2,
													  TLorentzVector& outLV_bu3) const {

	IMMMA_DecayResult4 result;

	//construct resonance2 = bu2 + bu3
	TLorentzVector resonance2_lab;

	//both bu2 and bu3 detected
	if(!fbu2.isMissing && !fbu3.isMissing){
		outLV_bu2 = BuildLab4Vector(fbu2);
		outLV_bu3 = BuildLab4Vector(fbu3);
		resonance2_lab = outLV_bu2 + outLV_bu3;
	}

	//if bu2 missing but bu3 detected
	else if(fbu2.isMissing && !fbu3.isMissing){
		outLV_bu3 = BuildLab4Vector(fbu3);

		//solve for bu2 from recoil and bu1
		if(!fbu1.isMissing){
			outLV_bu1 = BuildLab4Vector(fbu1);
			TVector3 p_res2 = recoil.Vect() - outLV_bu1.Vect();
			double E_res2 = recoil.E() - outLV_bu1.E();
			resonance2_lab.SetPxPyPzE(p_res2.X(), p_res2.Y(), p_res2.Z(), E_res2);

			//now get bu2 from resonance2 - bu3
			TVector3 p_bu2 = resonance2_lab.Vect() - outLV_bu3.Vect();
			double E_bu2 = resonance2_lab.E() - outLV_bu3.E();
			outLV_bu2.SetPxPyPzE(p_bu2.X(), p_bu2.Y(), p_bu2.Z(), E_bu2);
		} else {
			result.valid = false;
			return result;
		}
	}

	//if bu3 missing but bu2 detected
	else if(!fbu2.isMissing && fbu3.isMissing){
		outLV_bu2 = BuildLab4Vector(fbu2);

		//solve for bu3 from recoil and bu2
		if(!fbu1.isMissing){
			outLV_bu1 = BuildLab4Vector(fbu1);
			TVector3 p_res2 = recoil.Vect() - outLV_bu1.Vect();
			double E_res2 = recoil.E() - outLV_bu1.E();
			resonance2_lab.SetPxPyPzE(p_res2.X(), p_res2.Y(), p_res2.Z(), E_res2);

			//now get bu3 from resonance2 - bu2
			TVector3 p_bu3 = resonance2_lab.Vect() - outLV_bu2.Vect();
			double E_bu3 = resonance2_lab.E() - outLV_bu2.E();
			outLV_bu3.SetPxPyPzE(p_bu3.X(), p_bu3.Y(), p_bu3.Z(), E_bu3);
		} else {
			result.valid = false;
			return result;
		}
	}

	//both bu2 and bu3 missing
	else {
		result.valid = false;
		return result;
	}


	//switch to second decay resonance2->bu2+bu3
	// result.decay_resonance2 = SolveTwoBodyDecay(resonance2_lab, fbu2, fbu3, outLV_bu2, outLV_bu3);
	IMMMA_Fragment fbu2_analysis = fbu2;
	IMMMA_Fragment fbu3_analysis = fbu3;

	fbu2_analysis.isMissing = true;
	fbu3_analysis.isMissing = true;

	result.decay_resonance2 = SolveTwoBodyDecay(resonance2_lab, fbu2_analysis, fbu3_analysis, outLV_bu2, outLV_bu3);

	double res2_mass = resonance2_lab.M();
	double res2_KE = resonance2_lab.E() - res2_mass;
	double res2_theta = RADDEG * resonance2_lab.Theta();
	double res2_phi = RADDEG * resonance2_lab.Phi();
	if(res2_phi < 0) res2_phi += 360.;

	IMMMA_Fragment f_res2{res2_mass, res2_KE, res2_theta, res2_phi, false};



	//first decay
	IMMMA_Fragment fbu1_analysis = fbu1;
	IMMMA_Fragment f_res2_analysis = f_res2;
	TLorentzVector outLV_res2_dummy;
	fbu1_analysis.isMissing = true;
	f_res2_analysis.isMissing = true;

	result.decay_recoil = SolveTwoBodyDecay(recoil,
                      						fbu1_analysis,
                      						f_res2_analysis,
                      						outLV_bu1,
                      						outLV_res2_dummy);

	result.valid = result.decay_recoil.valid && result.decay_resonance2.valid;

	return result;
}

CaseResult4 IMMMA_Tool_4::ConvertDecayResult4(const IMMMA_DecayResult4& d,
											  const TLorentzVector& bu1_lab, double m1,
											  const TLorentzVector& bu2_lab, double m2,
											  const TLorentzVector& bu3_lab, double m3) const {

	CaseResult4 r;

	r.valid = d.valid;

	r.boostvector_recoil = d.decay_recoil.boost;

	r.Vcm_bu1_in_recoilCM = d.decay_recoil.Vcm1;
	r.KEcm_bu1_in_recoilCM = d.decay_recoil.KEcm1;
	r.ThetaCM_bu1_in_recoilCM = d.decay_recoil.ThetaCM1;
	r.PhiCM_bu1_in_recoilCM = d.decay_recoil.PhiCM1;

	r.Vcm_res2_in_recoilCM = d.decay_recoil.Vcm2;
	r.KEcm_res2_in_recoilCM = d.decay_recoil.KEcm2;
	r.ThetaCM_res2_in_recoilCM = d.decay_recoil.ThetaCM2;
	r.PhiCM_res2_in_recoilCM = d.decay_recoil.PhiCM2;

	r.Ecm_recoilCM = d.decay_recoil.Ecm;

	r.bu1_LabAngleWRTVCM_recoilCM = d.decay_recoil.LabAngle1WRTBoost;
	r.res2_LabAngleWRTVCM_recoilCM = d.decay_recoil.LabAngle2WRTBoost;

	r.boostvector_resonance2 = d.decay_resonance2.boost;

	r.Vcm_bu2_in_res2CM = d.decay_resonance2.Vcm1;
	r.KEcm_bu2_in_res2CM = d.decay_resonance2.KEcm1;
	r.ThetaCM_bu2_in_res2CM = d.decay_resonance2.ThetaCM1;
	r.PhiCM_bu2_in_res2CM = d.decay_resonance2.PhiCM1;

	r.Vcm_bu3_in_res2CM = d.decay_resonance2.Vcm2;
	r.KEcm_bu3_in_res2CM = d.decay_resonance2.KEcm2;
	r.ThetaCM_bu3_in_res2CM = d.decay_resonance2.ThetaCM2;
	r.PhiCM_bu3_in_res2CM = d.decay_resonance2.PhiCM2;

	r.Ecm_res2CM = d.decay_resonance2.Ecm;

	r.bu2_LabAngleWRTVCM_res2CM = d.decay_resonance2.LabAngle1WRTBoost;
	r.bu3_LabAngleWRTVCM_res2CM = d.decay_resonance2.LabAngle2WRTBoost;

	r.ELab_bu1 = bu1_lab.E() - m1;
	r.ELab_bu2 = bu2_lab.E() - m2;
	r.ELab_bu3 = bu3_lab.E() - m3;

	r.ThetaLab_bu1 = RADDEG * bu1_lab.Theta();
	r.PhiLab_bu1 = RADDEG * bu1_lab.Phi();
	if(r.PhiLab_bu1 < 0) r.PhiLab_bu1 += 360.;

	r.ThetaLab_bu2 = RADDEG * bu2_lab.Theta();
	r.PhiLab_bu2 = RADDEG * bu2_lab.Phi();
	if(r.PhiLab_bu2 < 0) r.PhiLab_bu2 += 360.;

	r.ThetaLab_bu3 = RADDEG * bu3_lab.Theta();
	r.PhiLab_bu3 = RADDEG * bu3_lab.Phi();
	if(r.PhiLab_bu3 < 0) r.PhiLab_bu3 += 360.;

	r.bu1InvMass = std::sqrt(std::abs(d.decay_recoil.invMassSquared1));
	r.res2InvMass = std::sqrt(std::abs(d.decay_recoil.invMassSquared2));
	r.recInvMass = std::sqrt(std::abs(d.decay_recoil.invMassSquaredParent));

	r.bu2InvMass = std::sqrt(std::abs(d.decay_resonance2.invMassSquared1));
	r.bu3InvMass = std::sqrt(std::abs(d.decay_resonance2.invMassSquared2));

	r.ejInvMass = 0.; // update this after returning otherwise it wont be passed to plot4mc

	return r;
}


std::array<CaseResult4, 6> IMMMA_Tool_4::AnalyzeEventIMM(double ejectileE, double ejectileTheta, double ejectilePhi,
														 double bu1E, double bu1Theta, double bu1Phi,
														 double bu2E, double bu2Theta, double bu2Phi,
														 double bu3E, double bu3Theta, double bu3Phi) const {


	if(breakups.size() < 3){
		throw std::runtime_error("IMMMA_Tool_4 requires 3 breakup nuclei");
	}

	//build recoil from all three detected break up particles
	TLorentzVector bu1temp = BuildLab4Vector(MakeFragment(breakups[0].massMeV, bu1E, bu1Theta, bu1Phi, false));
	TLorentzVector bu2temp = BuildLab4Vector(MakeFragment(breakups[1].massMeV, bu2E, bu2Theta, bu2Phi, false));
	TLorentzVector bu3temp = BuildLab4Vector(MakeFragment(breakups[2].massMeV, bu3E, bu3Theta, bu3Phi, false));

	TLorentzVector recoilLV = bu1temp + bu2temp + bu3temp;

	std::array<CaseResult4, 6> results;

	//6 permutations hard coded:							  indices:
	//perm 0: bu1=particle2, bu2=particle3, bu3=particle4 --> (0,1,2)
	//perm 1: bu1=particle2, bu2=particle4, bu3=particle3 --> (0,2,1)
	//perm 2: bu1=particle3, bu2=particle4, bu3=particle4 --> (1,0,2)
	//perm 3: bu1=particle3, bu2=particle4, bu3=particle2 --> (1,2,0)
	//perm 4: bu1=particle4, bu2=particle2, bu3=particle3 --> (2,0,1)
	//perm 5: bu1=particle4, bu2=particle3, bu3=particle2 --> (2,1,0)

	std::vector<std::array<int,3>> permutations = {
		{0,1,2},//perm 0
		{0,2,1},//perm 1
		{1,0,2},//perm 2
		{1,2,0},//perm 3
		{2,0,1},//perm 4
		{2,1,0} //perm 5
	};

	std::vector<std::array<double, 3>> detectedE = {
		{bu1E, bu2E, bu3E}
	};

	std::vector<std::array<double, 3>> detectedTheta = {
		{bu1Theta, bu2Theta, bu3Theta}
	};

	std::vector<std::array<double,3>> detectedPhi = {
		{bu1Phi, bu2Phi, bu3Phi}
	};


	for(size_t i = 0; i<6; i++){
		const auto& perm = permutations[i];

		IMMMA_Fragment fbu1 = MakeFragment(breakups[0],
											detectedE[0][perm[0]],
											detectedTheta[0][perm[0]],
											detectedPhi[0][perm[0]],
											false);

		IMMMA_Fragment fbu2 = MakeFragment(breakups[1],
											detectedE[0][perm[1]],
											detectedTheta[0][perm[1]],
											detectedPhi[0][perm[1]],
											false);

		IMMMA_Fragment fbu3 = MakeFragment(breakups[2],
											detectedE[0][perm[2]],
											detectedTheta[0][perm[2]],
											detectedPhi[0][perm[2]],
											false);

		TLorentzVector outLV_bu1, outLV_bu2, outLV_bu3;
		IMMMA_DecayResult4 decay = SolveSequentialDecay(recoilLV, fbu1, fbu2, fbu3, outLV_bu1, outLV_bu2, outLV_bu3);

		results[i] = ConvertDecayResult4(decay,
										 outLV_bu1, breakups[0].massMeV,
										 outLV_bu2, breakups[1].massMeV,
										 outLV_bu3, breakups[2].massMeV);

	}

	return results;

}


std::array<CaseResult4, 6> IMMMA_Tool_4::AnalyzeEventMMM(double ejectileE, double ejectileTheta, double ejectilePhi,
														 double detected1E, double detected1Theta, double detected1Phi,
														 double detected2E, double detected2Theta, double detected2Phi) const {


	if(breakups.size() < 3){
		throw std::runtime_error("IMMMA_Tool_4 requries 3 breakup nuclei");
	}

	TLorentzVector recoilLV = BuildBeamLV() + BuildTargetLV() - BuildEjectileLV(ejectileE, ejectileTheta, ejectilePhi);

	std::array<CaseResult4, 6> results;

	// There are 3 choices for which particle is missing * 2 ways to assign detected particles = 6 cases
	// Recall: resonance1/recoil -> breakup1 + resonance2, resonance2 -> breakup2 + breakup3 --> only can detect breakup1,2,3 and ejectile, resonance1/2 not directly measured
	// Hypothesis scheme (det1 and det2 are the two detected breakup particles):
	// hypothesis 0: det1=bu1, det2=bu2, missing=bu3
	// hypothesis 1: det1=bu2, det2=bu1, missing=bu3
	// hypothesis 2: det1=bu1, det2=bu3, missing=bu2
	// hypothesis 3: det1=bu3, det2=bu1, missing=bu2
	// hypothesis 4: det1=bu2, det2=bu3, missing=bu1
	// hypothesis 5: det1=bu3, det2=bu2, missing=bu1


	//-----------------------
	//--	hypothesis 0:  --
	//--	 bu1 && bu2    --
	//--	  detected,    --
	//--	 bu3 missed    --
	//-----------------------
	{
		IMMMA_Fragment fbu1 = MakeFragment(breakups[0], detected1E, detected1Theta, detected1Phi, false);
		IMMMA_Fragment fbu2 = MakeFragment(breakups[1], detected2E, detected2Theta, detected2Phi, false);
		IMMMA_Fragment fbu3 = MakeFragment(breakups[2], 0, 0, 0, true);
		TLorentzVector outLVbu1, outLVbu2, outLVbu3;
		IMMMA_DecayResult4 decay = SolveSequentialDecay(recoilLV, fbu1, fbu2, fbu3, outLVbu1, outLVbu2, outLVbu3);
		results[0] = ConvertDecayResult4(decay, outLVbu1, breakups[0].massMeV, outLVbu2, breakups[1].massMeV, outLVbu3, breakups[2].massMeV);
	}

	//-----------------------
	//--	hypothesis 1:  --
	//--	 bu2 && bu1    --
	//--	  detected,    --
	//--	 bu3 missed    --
	//-----------------------
	{
		IMMMA_Fragment fbu1 = MakeFragment(breakups[0], detected2E, detected2Theta, detected2Phi, false);
		IMMMA_Fragment fbu2 = MakeFragment(breakups[1], detected1E, detected1Theta, detected1Phi, false);
		IMMMA_Fragment fbu3 = MakeFragment(breakups[2], 0, 0, 0, true);
		TLorentzVector outLVbu1, outLVbu2, outLVbu3;
		IMMMA_DecayResult4 decay = SolveSequentialDecay(recoilLV, fbu1, fbu2, fbu3, outLVbu1, outLVbu2, outLVbu3);
		results[1] = ConvertDecayResult4(decay, outLVbu1, breakups[0].massMeV, outLVbu2, breakups[1].massMeV, outLVbu3, breakups[2].massMeV);
	}

	//-----------------------
	//--	hypothesis 2:  --
	//--	 bu1 && bu3    --
	//--	  detected,    --
	//--	 bu2 missed    --
	//-----------------------
	{
		IMMMA_Fragment fbu1 = MakeFragment(breakups[0], detected1E, detected1Theta, detected1Phi, false);
		IMMMA_Fragment fbu2 = MakeFragment(breakups[1], 0, 0, 0, true);
		IMMMA_Fragment fbu3 = MakeFragment(breakups[2], detected2E, detected2Theta, detected2Phi, false);
		TLorentzVector outLVbu1, outLVbu2, outLVbu3;
		IMMMA_DecayResult4 decay = SolveSequentialDecay(recoilLV, fbu1, fbu2, fbu3, outLVbu1, outLVbu2, outLVbu3);
		results[2] = ConvertDecayResult4(decay, outLVbu1, breakups[0].massMeV, outLVbu2, breakups[1].massMeV, outLVbu3, breakups[2].massMeV);
	}

	//-----------------------
	//--	hypothesis 3:  --
	//--	 bu3 && bu1    --
	//--	  detected,    --
	//--	 bu2 missed    --
	//-----------------------
	{
		IMMMA_Fragment fbu1 = MakeFragment(breakups[0], detected2E, detected2Theta, detected2Phi, false);
		IMMMA_Fragment fbu2 = MakeFragment(breakups[1], 0, 0, 0, true);
		IMMMA_Fragment fbu3 = MakeFragment(breakups[2], detected1E, detected1Theta, detected1Phi, false);
		TLorentzVector outLVbu1, outLVbu2, outLVbu3;
		IMMMA_DecayResult4 decay = SolveSequentialDecay(recoilLV, fbu1, fbu2, fbu3, outLVbu1, outLVbu2, outLVbu3);
		results[3] = ConvertDecayResult4(decay, outLVbu1, breakups[0].massMeV, outLVbu2, breakups[1].massMeV, outLVbu3, breakups[2].massMeV);
	}

	//-----------------------
	//--	hypothesis 4:  --
	//--	 bu2 && bu3    --
	//--	  detected,    --
	//--	 bu1 missed    --
	//-----------------------
	{
		IMMMA_Fragment fbu1 = MakeFragment(breakups[0], 0, 0, 0, true);
		IMMMA_Fragment fbu2 = MakeFragment(breakups[1], detected1E, detected1Theta, detected1Phi, false);
		IMMMA_Fragment fbu3 = MakeFragment(breakups[2], detected2E, detected2Theta, detected2Phi, false);
		TLorentzVector outLVbu1, outLVbu2, outLVbu3;
		IMMMA_DecayResult4 decay = SolveSequentialDecay(recoilLV, fbu1, fbu2, fbu3, outLVbu1, outLVbu2, outLVbu3);
		results[4] = ConvertDecayResult4(decay, outLVbu1, breakups[0].massMeV, outLVbu2, breakups[1].massMeV, outLVbu3, breakups[2].massMeV);
	}
	//-----------------------
	//--	hypothesis 5:  --
	//--	 bu3 && bu2    --
	//--	  detected,    --
	//--	 bu1 missed    --
	//-----------------------
	{
		IMMMA_Fragment fbu1 = MakeFragment(breakups[0], 0, 0, 0, true);
		IMMMA_Fragment fbu2 = MakeFragment(breakups[1], detected2E, detected2Theta, detected2Phi, false);
		IMMMA_Fragment fbu3 = MakeFragment(breakups[2], detected1E, detected1Theta, detected1Phi, false);
		TLorentzVector outLVbu1, outLVbu2, outLVbu3;
		IMMMA_DecayResult4 decay = SolveSequentialDecay(recoilLV, fbu1, fbu2, fbu3, outLVbu1, outLVbu2, outLVbu3);
		results[5] = ConvertDecayResult4(decay, outLVbu1, breakups[0].massMeV, outLVbu2, breakups[1].massMeV, outLVbu3, breakups[2].massMeV);
	}

	return results;
}
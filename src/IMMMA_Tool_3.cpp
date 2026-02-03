#include "IMMMA_Tool_3.h"

IMMMA_Tool_3::IMMMA_Tool_3() = default;

//setters:
void IMMMA_Tool_3::SetBeamNucleus(int A, TString sym, double mass){
	beam = {A, sym, mass};
}

void IMMMA_Tool_3::SetTargetNucleus(int A, TString sym, double mass){
	target = {A, sym, mass};
}

void IMMMA_Tool_3::SetEjectileNucleus(int A, TString sym, double mass){
	ejectile = {A, sym, mass};
}

void IMMMA_Tool_3::SetRecoilNucleus(int A, TString sym, double mass){
	recoil = {A, sym, mass};
}

void IMMMA_Tool_3::SetBreakupNuclei(const std::vector<Nucleus>& b){
	breakups = b;
}

void IMMMA_Tool_3::SetBeamEnergyMeV(double energy){
	beam_energy = energy;
}

void IMMMA_Tool_3::SetRecoilExEMeV(double energy){
	recoil_ExE = energy;
}

//build lab-frame vectors:
TLorentzVector IMMMA_Tool_3::BuildBeamLV() const {
	// TLorentzVector v;
	// double p = std::sqrt(2.0 * beam.massMeV * beam_energy);
	// v.SetPxPyPzE(0, 0, p, beam.massMeV + beam_energy);
	// return v;
	IMMMA_Fragment f{beam.massMeV, beam_energy, 0, 0, false};
	return BuildLab4Vector(f);
}

TLorentzVector IMMMA_Tool_3::BuildTargetLV() const {
	// TLorentzVector v;
	// v.SetPxPyPzE(0, 0, 0, target.massMeV);
	// return v;
	IMMMA_Fragment f{target.massMeV, 0, 0, 0, false};
	return BuildLab4Vector(f);
}

TLorentzVector IMMMA_Tool_3::BuildEjectileLV(double E, double th, double ph) const {
	// TLorentzVector v;
	// double p = std::sqrt(2.0 * ejectile.massMeV * E);

	// v.SetPxPyPzE(
	// 				p*std::sin(DEGRAD * th) * std::cos(DEGRAD * ph),
	// 				p*std::sin(DEGRAD * th) * std::sin(DEGRAD * ph),
	// 				p*std::cos(DEGRAD * th),
	// 				ejectile.massMeV + E
	// 			);
	// return v;

	IMMMA_Fragment f{ejectile.massMeV, E, th, ph, false};
	return BuildLab4Vector(f);
}

CaseResult IMMMA_Tool_3::ConvertDecayResult(const IMMMA_DecayResult& d,
											const TLorentzVector& p1_lab, double m1,
											const TLorentzVector& p2_lab, double m2
										   ) const {


	CaseResult r;

	r.boostvector = d.boost;

	r.Vcm1 = d.Vcm1; r.KEcm1 = d.KEcm1; r.ThetaCM1 = d.ThetaCM1; r.PhiCM1 = d.PhiCM1;
	r.Vcm2 = d.Vcm2; r.KEcm2 = d.KEcm2; r.ThetaCM2 = d.ThetaCM2; r.PhiCM2 = d.PhiCM2;

	r.Ecm = d.Ecm;

	r.breakup1_LabAngleWRTVCM = d.LabAngle1WRTBoost;
	r.breakup2_LabAngleWRTVCM = d.LabAngle2WRTBoost;

	r.ELab1 = p1_lab.Energy();
	r.ELab2 = p2_lab.Energy();

	r.ThetaLab1 = RADDEG * p1_lab.Theta();
	r.PhiLab1 = RADDEG * p1_lab.Phi();
	if(r.PhiLab1 < 0) r.PhiLab1 += 360.;

	r.ThetaLab2 = RADDEG * p2_lab.Theta();
	r.PhiLab2 = RADDEG * p2_lab.Phi();
	if(r.PhiLab2 < 0) r.PhiLab2 += 360.;

	r.ejInvMass = std::sqrt(d.invMassSquaredEj);
	r.bu1InvMass = std::sqrt(d.invMassSquared1);
	r.bu2InvMass = std::sqrt(d.invMassSquared2);
	r.recInvMass = std::sqrt(d.invMassSquaredParent);

	return r;
}


// IMM Analysis:
std::pair<CaseResult, CaseResult> IMMMA_Tool_3::AnalyzeEventIMM(
		double ejectileE, double ejectileTheta, double ejectilePhi,
		double detected1E, double detected1Theta, double detected1Phi,
		double detected2E, double detected2Theta, double detected2Phi) const
{

	if(breakups.size() <= 2){
		std::cout << breakups.size() << std::endl;
		throw std::runtime_error("IMMMA_Tool_3 requires 2 breakup nuclei!");
	}

	//TLorentzVector recoilLV = BuildBeamLV() + BuildTargetLV() - BuildEjectileLV(ejectileE, ejectileTheta, ejectilePhi);
	TLorentzVector recoilLV = BuildLab4Vector(MakeFragment(breakups[0].massMeV, detected1E, detected1Theta, detected1Phi, false)) + BuildLab4Vector(MakeFragment(breakups[1].massMeV, detected2E, detected2Theta, detected2Phi, false));

	//hypothesis A:
	IMMMA_Fragment f1A = MakeFragment(breakups[0], detected1E, detected1Theta, detected1Phi);
	IMMMA_Fragment f2A = MakeFragment(breakups[1], detected2E, detected2Theta, detected2Phi);

	TLorentzVector outLV1A, outLV2A;
	IMMMA_DecayResult dA = SolveTwoBodyDecay(recoilLV, f1A, f2A, outLV1A, outLV2A);
	CaseResult resultA = ConvertDecayResult(dA, outLV1A, breakups[0].massMeV, outLV2A, breakups[1].massMeV);

	//hypothesis B:
	IMMMA_Fragment f1B = MakeFragment(breakups[1], detected1E, detected1Theta, detected1Phi);
	IMMMA_Fragment f2B = MakeFragment(breakups[0], detected2E, detected2Theta, detected2Phi);
	TLorentzVector outLV1B, outLV2B;
	IMMMA_DecayResult dB = SolveTwoBodyDecay(recoilLV, f1B, f2B, outLV1B, outLV2B);
	CaseResult resultB = ConvertDecayResult(dB, outLV1B, breakups[1].massMeV, outLV2B, breakups[0].massMeV);

	return{resultA, resultB};

}

// MMM Analysis:
std::pair<CaseResult, CaseResult> IMMMA_Tool_3::AnalyzeEventMMM(
													double ejectileE, double ejectileTheta, double ejectilePhi,
													double detectedE, double detectedTheta, double detectedPhi) const
{

	if(breakups.size() <= 2){
		std::cout << breakups.size() << std::endl;
		throw std::runtime_error("IMMMA_Tool_3 requires 2 breakup nuclei!");
	}


	TLorentzVector recoilLV = BuildBeamLV() + BuildTargetLV() - BuildEjectileLV(ejectileE, ejectileTheta, ejectilePhi);

	// hypothesis A:
	IMMMA_Fragment obsA = MakeFragment(breakups[0], detectedE, detectedTheta, detectedPhi);
	IMMMA_Fragment missA = MakeFragment(breakups[1], 0, 0, 0, true);
	TLorentzVector outLV1A, outLV2A;
	IMMMA_DecayResult decayA = SolveTwoBodyDecay(recoilLV, obsA, missA, outLV1A, outLV2A);
	CaseResult resultA = ConvertDecayResult(decayA, outLV1A, breakups[0].massMeV, outLV2A, breakups[1].massMeV);

	// hypothesis  B:
	IMMMA_Fragment obsB = MakeFragment(breakups[1], detectedE, detectedTheta, detectedPhi);
	IMMMA_Fragment missB = MakeFragment(breakups[0], 0, 0, 0, true);
	TLorentzVector outLV1B, outLV2B;
	IMMMA_DecayResult decayB = SolveTwoBodyDecay(recoilLV, obsB, missB, outLV1B, outLV2B);
	CaseResult resultB = ConvertDecayResult(decayB, outLV1B, breakups[1].massMeV, outLV2B, breakups[0].massMeV);

	return {resultA, resultB};
}
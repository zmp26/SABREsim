#ifndef IMMMA_TOOL_3_H
#define IMMMA_TOOL_3_H

#include "IMMMA_Tool_Base.h"
#include "structs.h"

#include "TString.h"
#include "TLorentzVector.h"
#include <utility>


class IMMMA_Tool_3 : protected IMMMA_Tool_Base {
public:
	IMMMA_Tool_3();

	//setters
	void SetBeamNucleus(int A, TString sym, double mass);
	void SetTargetNucleus(int A, TString sym, double mass);
	void SetEjectileNucleus(int A, TString sym, double mass);
	void SetRecoilNucleus(int A, TString sym, double mass);
	void SetBreakup1Nucleus(int A, TString sym, double mass);
	void SetBreakup2Nucleus(int A, TString sym, double mass);

	void SetBeamEnergyMeV(double energy);
	void SetRecoilExEMeV(double energy);


	//analysis
	std::pair<CaseResult, CaseResult> AnalyzeEventIMM(
														double ejectileE, double ejectileTheta, double ejectilePhi,
														double detected1E, double detected1Theta, double detected1Phi,
														double detected2E, double detected2Theta, double detected2Phi
													 ) const;


	std::pair<CaseResult, CaseResult> AnalyzeEventMMM(
														double ejectileE, double ejectileTheta, double ejectilePhi,
														double detected1E, double detected1Theta, double detected1Phi
													 ) const;


private:
	Nucleus beam;
	Nucleus target;
	Nucleus ejectile;
	Nucleus recoil;
	Nucleus breakup1;
	Nucleus breakup2;

	double beam_energy = 0;
	double recoil_ExE = 0;

	CaseResult ConvertDecayResult(const IMMMA_DecayResult& d, const TLorentzVector& p1_lab, const TLorentzVector& p2_lab) const;

	TLorentzVector BuildBeamLV() const;
	TLorentzVector BuildTargetLV() const;
	TLorentzVector BuildEjectileLV(double E, double th, double ph) const;
};

#endif//IMMMA_TOOL_3_H
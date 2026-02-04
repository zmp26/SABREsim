#ifndef IMMMA_TOOL_4_H
#define IMMMA_TOOL_4_H

#include "IMMMA_Tool_Base.h"
#include "structs.h"
#include <vector>
#include <array>

// struct CaseResult4 {

// 	//boost vectors
// 	TVector3 boostvector_recoil;			//VCM_recoil = VCM_resonance1		resonance1 -> resonance2 + breakup1
// 	TVector3 boostvector_resonance2;		//VCM_resonance2					resonance2 -> breakup2 + breakup3


// 	//--------------------------------------
// 	//-		recoil/resonance1 CM frame	   -
// 	//--------------------------------------

// 	//breakup1 in recoil/resonance1 CM frame
// 	double Vcm_bu1_in_recoilCM;
// 	double KEcm_bu1_in_recoilCM;
// 	double ThetaCM_bu1_in_recoilCM;
// 	double PhiCM_bu1_in_recoilCM;

// 	//resonance2 in recoil/resonance1 CM frame (this is not measured)
// 	double Vcm_res2_in_recoilCM;
// 	double KEcm_res2_in_recoilCM;
// 	double ThetaCM_res2_in_recoilCM;
// 	double PhiCM_res2_in_recoilCM;

// 	//total energy in recoil/resonance1 CM
// 	double Ecm_recoilCM;



// 	//--------------------------------------
// 	//-		recoil/resonance1 CM frame	   -
// 	//--------------------------------------

// 	//breakup2 in resonance2 CM frame
// 	double Vcm_bu2_in_res2CM;
// 	double KEcm_bu2_in_res2CM;
// 	double ThetaCM_bu2_in_res2CM;
// 	double PhiCM_bu2_in_res2CM;

// 	//breakup3 in resonance2 CM frame
// 	double Vcm_bu3_in_res2CM;
// 	double KEcm_bu3_in_res2CM;
// 	double ThetaCM_bu3_in_res2CM;
// 	double PhiCM_bu3_in_res2CM;

// 	//total energy in resonance2 CM
// 	double Ecm_res2CM;

// 	//--------------------------------------
// 	//-		Lab frame angles WRT VCM	   -
// 	//--------------------------------------
// 	double bu1_LabAngleWRTVCM_recoilCM;
// 	double res2_LabAngleWRTVCM_recoilCM;
// 	double breakup2_LabAngleWRTVCM_res2CM;
// 	double breakup2_LabAngleWRTVCM_res2CM;

// 	//--------------------------------------
// 	//-			Lab frame energies 		   -
// 	//--------------------------------------
// 	double ELab_bu1;
// 	double ELab_bu2;
// 	double ELab_bu3;

// 	//--------------------------------------
// 	//-			Lab frame angles 		   -
// 	//--------------------------------------
// 	double ThetaLab_bu1;
// 	double PhiLab_bu1;
// 	double ThetaLab_bu2;
// 	double PhiLab_bu2;
// 	double ThetaLab_bu3;
// 	double PhiLab_bu3;


// 	//--------------------------------------
// 	//-			Invariant Masses 		   -
// 	//--------------------------------------
// 	double ejInvMass;
// 	double bu1InvMass;
// 	double bu2InvMass;
// 	double bu3InvMass;
// 	double res2InvMass;
// 	double recInvMass;


// 	//--------------------------------------
// 	//-			  Validity Flag		   	   -
// 	//--------------------------------------
// 	bool valid = true;
// };

class IMMMA_Tool_4 : public IMMMA_Tool_Base {
private:
	Nucleus beam;
	Nucleus target;
	Nucleus ejectile;
	Nucleus recoil;
	std::vector<Nucleus> breakups; //should be size() = 3 (breakup1, breakup2, breakup3)

	double beam_energy = 0.0;
	double recoil_ExE = 0.0;

	TLorentzVector BuildBeamLV() const;
	TLorentzVector BuildTargetLV() const;
	TLorentzVector BuildEjectileLV(double E, double th, double ph) const;

	CaseResult4 ConvertDecayResult4(const IMMMA_DecayResult4& d,
								   const TLorentzVector& bu1_lab, double m1,
								   const TLorentzVector& bu2_lab, double m2,
								   const TLorentzVector& bu3_lab, double m3) const;

	IMMMA_DecayResult4 SolveSequentialDecay(const TLorentzVector& recoil,
											const IMMMA_Fragment& fbu1,
											const IMMMA_Fragment& fbu2,
											const IMMMA_Fragment& fbu3,
											TLorentzVector& outLV_bu1,
											TLorentzVector& outLV_bu2,
											TLorentzVector& outLV_bu3) const;

public:
	IMMMA_Tool_4();

	//setters
	void SetBeamNucleus(int A, TString sym, double mass);
	void SetTargetNucleus(int A, TString sym, double mass);
	void SetEjectileNucleus(int A, TString sym, double mass);
	void SetRecoilNucleus(int A, TString sym, double mass);
	void SetBreakupNuclei(const std::vector<Nucleus>& b); //expects 3 nuclei!
	void SetBeamEnergyMeV(double energy);
	void SetRecoilExEMeV(double energy);


	//analysis methods
	std::array<CaseResult4, 6> AnalyzeEventIMM(double ejectileE, double ejectileTheta, double ejectilePhi,
											   double bu1E, double bu1Theta, double bu1Phi,
											   double bu2E, double bu2Theta, double bu2Phi,
											   double bu3E, double bu3Theta, double bu3Phi) const;


	std::array<CaseResult4, 6> AnalyzeEventMMM(double ejectileE, double ejectileTheta, double ejectilePhi,
											   double detected1E, double detected1Theta, double detected1Phi,
											   double detected2E, double detected2Theta, double detected2Phi) const;

};

#endif//IMMMA_TOOL_4_H
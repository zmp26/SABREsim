#ifndef STRUCTS_H
#define STRUCTS_H

#include <cmath>
#include "TString.h"
#include "TVector3.h"

static const int numwedges = 8;
static const int numrings = 16;

static constexpr double DEGRAD = M_PI/180.;
static constexpr double RADDEG = 180./M_PI;

static const std::pair<int,int> offsets[] = {
	{112,40},	//detector0 {ringOffset,wedgeOffset}
	{96,32},	//detector1 {ringOffset,wedgeOffset}
	{80,16},	//detector2 {ringOffset,wedgeOffset}
	{64,24},	//detector3 {ringOffset,wedgeOffset}
	{48,0}		//detector4 {ringOffset,wedgeOffset}
};

struct PHYSDATA { double e, theta, phi; };
struct SABREDATA { int detectorIndex=-666, particleIndex=-666; double theta, phi, ringEnergy, wedgeEnergy, localx, localy; int ring, wedge; };

struct Nucleus {
	int A = -666;
	TString sym = "";
	double massMeV = -666;

	Nucleus() = default;
	Nucleus(int a, const TString& s, double m) : A(a), sym(s), massMeV(m) {}

	Nucleus(const Nucleus&) = default;

	Nucleus(Nucleus&&) = default;

	Nucleus& operator=(const Nucleus&) = default;

	bool Filled() const {
		return (A != -666 && sym != "" && massMeV != -666);
	}

	TString ToString() const {
		return Form("%d%s", A, sym.Data());
	}

	void SetA(int newA){
		A = newA;
	}

	void SetSym(TString newsym){
		sym = newsym;
	}

	void SetMassMeV(double newmassMeV){
		massMeV = newmassMeV;
	}

	void SetAll(int newA, TString newsym, double newmassMeV){
		SetA(newA);
		SetSym(newsym);
		SetMassMeV(newmassMeV);
	}
};

struct Reaction2{
	Nucleus beam;
	Nucleus target;
	Nucleus ejectile;
	Nucleus recoil;

	double beamEnergy=0.;
	double recoilExE=0.;

	TString ToString(){
		TString retval = Form("%d%s(%d%s,%d%s)%d%s at E=%f to ExE=%f",target.A,target.sym.Data(),beam.A,beam.sym.Data(),ejectile.A,ejectile.sym.Data(),recoil.A,recoil.sym.Data(),beamEnergy,recoilExE);
		return retval;
	};
};

struct Reaction3{
	Nucleus beam;
	Nucleus target;
	Nucleus ejectile;
	Nucleus recoil;
	Nucleus breakup1;
	Nucleus breakup2;

	double beamEnergy=0.;
	double recoilExE=0.;

	TString ToString(){
		TString retval = Form("%d%s(%d%s,%d%s)%d%s at E=%f to ExE=%f",target.A,target.sym.Data(),beam.A,beam.sym.Data(),ejectile.A,ejectile.sym.Data(),recoil.A,recoil.sym.Data(),beamEnergy,recoilExE);
		return retval;
	};
};

struct CaseResult {

	TVector3 boostvector;

	double Vcm1 = 0;
	double KEcm1 = 0;
	double ThetaCM1 = 0;
	double PhiCM1 = 0;
	double ELab1 = 0;
	double ThetaLab1 = 0;
	double PhiLab1 = 0;
	double breakup1_LabAngleWRTVCM = 0;

	double Vcm2 = 0;
	double KEcm2 = 0;
	double ThetaCM2 = 0;
	double PhiCM2 = 0;
	double ELab2 = 0;
	double ThetaLab2 = 0;
	double PhiLab2 = 0;
	double breakup2_LabAngleWRTVCM = 0;

	double recInvMass = 0;
	double ejInvMass = 0;
	double bu1InvMass = 0;
	double bu2InvMass = 0;

	double Ecm = 0;
	double recoilExE = 0;

};

struct CaseResult4 {

	//boost vectors
	TVector3 boostvector_recoil;			//VCM_recoil = VCM_resonance1		resonance1 -> resonance2 + breakup1
	TVector3 boostvector_resonance2;		//VCM_resonance2					resonance2 -> breakup2 + breakup3


	//--------------------------------------
	//-		recoil/resonance1 CM frame	   -
	//--------------------------------------

	//breakup1 in recoil/resonance1 CM frame
	double Vcm_bu1_in_recoilCM;
	double KEcm_bu1_in_recoilCM;
	double ThetaCM_bu1_in_recoilCM;
	double PhiCM_bu1_in_recoilCM;

	//resonance2 in recoil/resonance1 CM frame (this is not measured)
	double Vcm_res2_in_recoilCM;
	double KEcm_res2_in_recoilCM;
	double ThetaCM_res2_in_recoilCM;
	double PhiCM_res2_in_recoilCM;

	//total energy in recoil/resonance1 CM
	double Ecm_recoilCM;



	//--------------------------------------
	//-		recoil/resonance1 CM frame	   -
	//--------------------------------------

	//breakup2 in resonance2 CM frame
	double Vcm_bu2_in_res2CM;
	double KEcm_bu2_in_res2CM;
	double ThetaCM_bu2_in_res2CM;
	double PhiCM_bu2_in_res2CM;

	//breakup3 in resonance2 CM frame
	double Vcm_bu3_in_res2CM;
	double KEcm_bu3_in_res2CM;
	double ThetaCM_bu3_in_res2CM;
	double PhiCM_bu3_in_res2CM;

	//total energy in resonance2 CM
	double Ecm_res2CM;

	//--------------------------------------
	//-		Lab frame angles WRT VCM	   -
	//--------------------------------------
	double bu1_LabAngleWRTVCM_recoilCM;
	double res2_LabAngleWRTVCM_recoilCM;
	double bu2_LabAngleWRTVCM_res2CM;
	double bu3_LabAngleWRTVCM_res2CM;

	//--------------------------------------
	//-			Lab frame energies 		   -
	//--------------------------------------
	double ELab_bu1;
	double ELab_bu2;
	double ELab_bu3;

	//--------------------------------------
	//-			Lab frame angles 		   -
	//--------------------------------------
	double ThetaLab_bu1;
	double PhiLab_bu1;
	double ThetaLab_bu2;
	double PhiLab_bu2;
	double ThetaLab_bu3;
	double PhiLab_bu3;


	//--------------------------------------
	//-			Invariant Masses 		   -
	//--------------------------------------
	double ejInvMass;
	double bu1InvMass;
	double bu2InvMass;
	double bu3InvMass;
	double res2InvMass;
	double recInvMass;


	//--------------------------------------
	//-			  Validity Flag		   	   -
	//--------------------------------------
	bool valid = true;
};

#endif//STRUCTS_H
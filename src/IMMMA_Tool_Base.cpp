#include "IMMMA_Tool_Base.h"

TLorentzVector IMMMA_Tool_Base::BuildLab4Vector(const IMMMA_Fragment& f) const{

	TLorentzVector v;

	if(f.isMissing){
		v.SetPxPyPzE(0., 0., 0., f.mass);
		return v;
	}

	double p = std::sqrt(2.0 * f.mass * f.Elab);

	v.SetPxPyPzE(
		p*std::sin(DEGRAD * f.theta) * std::cos(DEGRAD * f.phi),
		p*std::sin(DEGRAD * f.theta) * std::sin(DEGRAD * f.phi),
		p*std::cos(DEGRAD * f.theta),
		f.mass + f.Elab
	);

	return v;
}

IMMMA_Fragment IMMMA_Tool_Base::MakeFragment(const Nucleus& nuc, double E, double th, double ph, bool missing) const {
	return IMMMA_Fragment{nuc.massMeV, E, th, ph, missing};
}

IMMMA_DecayResult IMMMA_Tool_Base::SolveTwoBodyDecay(
	const TLorentzVector& parent,
	const IMMMA_Fragment& f1,
	const IMMMA_Fragment& f2,
	TLorentzVector& outLV1,
	TLorentzVector& outLV2
) const {

	IMMMA_DecayResult res;


	//compute fragment 4-vectors
	//case 1: both fragments detected
	if(!f1.isMissing && !f2.isMissing){
		outLV1 = BuildLab4Vector(f1);
		outLV2 = BuildLab4Vector(f2);
	}

	//case 2: f1 missing, f2 detected
	else if(f1.isMissing && !f2.isMissing){
		outLV2 = BuildLab4Vector(f2);
		outLV1 = parent - outLV2;
	}

	//case 3: f1 detected, f2 missing
	else if(!f1.isMissing && f2.isMissing){
		outLV1 = BuildLab4Vector(f1);
		outLV2 = parent - outLV1;
	}

	//case 4: neither detected
	else{
		outLV1.SetPxPyPzE(0,0,0,f1.mass);
		outLV2.SetPxPyPzE(0,0,0,f2.mass);
		res.valid = false;
		return res;
	}


	//compute composite lorentz vector and boost vector:
	TLorentzVector composite = outLV1 + outLV2;
	TVector3 boost = -composite.BoostVector();
	res.boost = boost;

	//compute lab angles WRT boost vector (VCM) in lab frame
	res.LabAngle1WRTBoost = RADDEG * outLV1.Vect().Angle(-boost);
	res.LabAngle2WRTBoost = RADDEG * outLV2.Vect().Angle(-boost);

	//compute CM frame kinematics
	TLorentzVector LV1_cm = outLV1;
	TLorentzVector LV2_cm = outLV2;
	LV1_cm.Boost(boost);
	LV2_cm.Boost(boost);

	//fragment 1 CM
	res.Vcm1 = (LV1_cm.Vect() * ( 1 / LV1_cm.E())).Mag();
	res.KEcm1 = 0.5 * f1.mass * res.Vcm1 * res.Vcm1;
	res.ThetaCM1 = RADDEG * LV1_cm.Theta();
	res.PhiCM1 = RADDEG * LV1_cm.Phi();
	if(res.PhiCM1 < 0) res.PhiCM1 += 360.;

	//fragment 2 CM
	res.Vcm2 = (LV2_cm.Vect() * ( 1 / LV2_cm.E())).Mag();
	res.KEcm2 = 0.5 * f2.mass * res.Vcm2 * res.Vcm2;
	res.ThetaCM2 = RADDEG * LV2_cm.Theta();
	res.PhiCM2 = RADDEG * LV2_cm.Phi();
	if(res.PhiCM2 < 0) res.PhiCM2 += 360.;

	res.Ecm = res.KEcm1 + res.KEcm2;

	//get squared invariant masses:
	res.invMassSquared1 = outLV1.M2();
	res.invMassSquared2 = outLV2.M2();
	res.invMassSquaredParent = composite.M2();

	return res;
}
#ifndef IMMMA_TOOL_BASE_H
#define IMMMA_TOOL_BASE_H

#include "structs.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <vector>
#include <cmath>

/*
 * IMMMA_Tool_Base
 *
 * Written by Zach Purcell zpurce2@lsu.edu
 *
 * ROOT-backed kinematics engine for reconstructing resonances through
 * invariant mass (IM) and missing mass (MM) analysis of multi-body breakups
 *
 * IMMMA_Tool_Base:
 *  - Handles arbitrary decay depth via repeated 2-body decays
 *  - Supports detected and missing fragments
 *
 *
 * IMMMA_Tool_Base contains NO experiment-specific logic!
 *
 */

struct IMMMA_Fragment {		// UNITS:
	double mass = 0.;		// MeV/c^2
	double Elab = 0.;		// MeV
	double theta = 0.;		// degrees
	double phi = 0.;		// degrees
	bool isMissing = false;
};

struct IMMMA_DecayResult {
	bool valid = true;//default - change to false if invalid for some reason (for example, if both fragments somehow marked as isMissing = true)

	TVector3 boost;

	//fragment 1 (passed in)
	double Vcm1 = 0.;
	double KEcm1 = 0.;
	double ThetaCM1 = 0.;
	double PhiCM1 = 0.;
	double LabAngle1WRTBoost = 0.;

	//fragment 2 (passed in)
	double Vcm2 = 0.;
	double KEcm2 = 0.;
	double ThetaCM2 = 0.;
	double PhiCM2 = 0.;
	double LabAngle2WRTBoost = 0.;

	//CM energies:
	double Ecm = 0.;

	//invariant masses squared:
	double invMassSquaredEj = 0.;
	double invMassSquared1 = 0.;
	double invMassSquared2 = 0.;
	double invMassSquaredParent = 0.;
};

struct IMMMA_EventResult{
	std::vector<IMMMA_DecayResult> decays;
	double recoilExE = 0.;
};

class IMMMA_Tool_Base{
public:
	IMMMA_Tool_Base() = default;
	virtual ~IMMMA_Tool_Base() = default;

protected:
	//core kinematics engine
	// IMMMA_DecayResult SolveTwoBodyDecay(
	// 	const TLorentzVector& parent,
	// 	const IMMMA_Fragment& f1,
	// 	const IMMMA_Fragment& f2,
	// 	TLorentzVector& outComposite
	// ) const;
	IMMMA_DecayResult SolveTwoBodyDecay(
		const TLorentzVector& parent,
		const IMMMA_Fragment& f1,
		const IMMMA_Fragment& f2,
		TLorentzVector& outLV1,
		TLorentzVector& outLV2
	) const;

	TLorentzVector BuildLab4Vector(const IMMMA_Fragment& f) const;

	IMMMA_Fragment MakeFragment(const Nucleus& nuc, double E, double th, double ph, bool missing = false) const;
	IMMMA_Fragment MakeFragment(double massMeV, double E, double th, double ph, bool missing = false) const;
};

#endif//IMMMA_TOOL_BASE_H
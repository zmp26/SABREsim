#ifndef PID_STRUCTS_H
#define PID_STRUCTS_H

#include <array>
#include "TLorentzVector.h"
#include "TString.h"

constexpr double DEGRAD = M_PI / 180.;
constexpr double RADDEG = 180. / M_PI;
constexpr double AMU_IN_MEV = 931.5;

// Shared N=3 Reaction Hypothesis Configuration
struct PIDHypothesis {
	TString name;
	double beamEnergyMeV{0.0};
	double mass_beam{0.0};
	double mass_target{0.0};
	double mass_ejectile{0.0};
	double final_masses[3]{0.0, 0.0, 0.0};     // Rest masses for species 0, 1, 2
	TString final_particles[3]{"", "", ""};    // Names/labels (e.g. "p", "alpha", "alpha")
};

// Shared N=2 Reaction Hypothesis Configuration
struct PIDHypothesis_N2{
	TString name;
	double beamEnergyMeV{0.};
	double mass_beam{0.0};
	double mass_target{0.0};
	double mass_ejectile{0.0};
	double final_masses[2]{0.,0.};
	TString final_particles[2]{"", ""};
};

// Result Container for Complete Kinematics (3 Hits Detected)
//TODO: Update PIDResult_Mult3 to PIDResult_N3_M3
struct PIDResult_Mult3 {
	int bestChi2Index{-1};
	double bestChi2{1e9};
	bool passesCut{false};

	std::array<int, 3> hit_indices{-1, -1, -1}; // Maps final_masses[0,1,2] species to detected hits H0, H1, H2

	// Momentum residual vector (P_recoil - sum(P_decay))
	double missing_px{0.0};
	double missing_py{0.0};
	double missing_pz{0.0};
	double missing_E{0.0};
	double missing_Pmag{0.0};

	double m2_aa{0.};
	double m2_pa1{0.};
	double m2_pa2{0.};
	double ExSPS{0.};

	std::array<double, 6> permChi2s;

	void Reset() {
		bestChi2Index = -1;
		bestChi2 = 1e9;
		passesCut = false;
		hit_indices = {-1, -1, -1};
		missing_px = 0.0;
		missing_py = 0.0;
		missing_pz = 0.0;
		missing_E = 0.0;
		missing_Pmag = 0.0;
		permChi2s.fill(1e9);
		m2_aa = 0.;
		m2_pa1 = 0.;
		m2_pa2 = 0.;
		ExSPS = 0.;
	}
};

// Result Container for Incomplete Kinematics (2 Hits Detected, 1 Missing)
//TODO: Update PIDResult_Mult2 to PIDResult_N3_M2
struct PIDResult_Mult2 {
	int bestChi2Index{-1};
	double bestChi2{1e9};
	bool passesCut{false};

	std::array<int, 2> hit_indices{-1, -1};  // Species indices (0, 1, or 2) assigned to detected hits H0, H1
	int missing_species_index{-1};          // Species index (0, 1, or 2) assigned to undetected particle

	// Physical reconstructed 4-momentum of the missing particle
	double missing_px{0.0};
	double missing_py{0.0};
	double missing_pz{0.0};
	double missing_Etot{0.0};
	double missing_Pmag{0.0};
	double missing_MassCalc{0.0};            // Reconstructed invariant mass: sqrt(P4_missing^2)

	std::array<double, 6> permChi2s;

	void Reset() {
		bestChi2Index = -1;
		bestChi2 = 1e9;
		passesCut = false;
		hit_indices = {-1, -1};
		missing_species_index = -1;
		missing_px = 0.0;
		missing_py = 0.0;
		missing_pz = 0.0;
		missing_Etot = 0.0;
		missing_Pmag = 0.0;
		missing_MassCalc = 0.0;
		permChi2s.fill(1e9);
	}
};

// result container for complete kinematics in N=2
struct PIDResult_N2_M2 {
	double bestChi2{1e9};
	int bestChi2Index{-1};

	std::array<double,2> permChi2s{{1e9,1e9}};
	std::array<int,2> hit_indices{{-1,-1}};

	double missing_px{0.};
	double missing_py{0.};
	double missing_pz{0.};
	double missing_E{0.};
	double missing_Pmag{0.};

	bool passesCut{false};

	void Reset(){
		bestChi2 = 1e9;
		bestChi2Index = -1;
		permChi2s = {{1e9, 1e9}};
		hit_indices = {{-1, -1}};

		missing_px = 0.0;
		missing_py = 0.0;
		missing_pz = 0.0;
		missing_E = 0.0;
		missing_Pmag = 0.0;

		passesCut = false;
	}

};

double Kallen(double a, double b, double c){
	//the Kallen function, or triangle function, often denoted lambda(a,b,c) = a^2 + b^2 + c^2 - 2ab - 2ac - 2bc
	return a*a + b*b + c*c - 2.*(a*b + a*c + b*c);
}

#endif // PID_STRUCTS_H
#ifndef INVMASSMULT2_H
#define INVMASSMULT2_H

#include <map>
#include <vector>
#include <array>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
//#include "permHisto_mult2.h"

//Hypothesis3 struct
struct Hypothesis3 {
	std::string name;

	double mass_target;
	double mass_beam;
	double mass_ejectile;
	double mass_recoil;
	double mass_intermediate;

	double masses[2];

};

enum gateIndices {
	NOCHECK,
	FIRSTVALID = NOCHECK,
	INTEXCHECK,
	INTVCMCHECK,
	FRAG1VCMCHECK,
	LASTVALID = FRAG
}


#endif//INVMASSMULT2_H
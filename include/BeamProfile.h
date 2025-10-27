#ifndef BEAM_PROFILE_H
#define BEAM_PROFILE_H

#include <utility>
#include <random>
#include "TString.h"

class BeamProfile{
public:
	virtual ~BeamProfile() = default;

	virtual std::pair<double, double> Sample() = 0;

	virtual TString ToString() const { return "BeamProfile"; };

	double GetParX() { return parX; };
	double GetParY() { return parY; };

protected:
	double parX=0.;
	double parY=0.;

};

#endif //BEAM_PROFILE_H
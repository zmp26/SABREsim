#ifndef FIXED_POINT_PROFILE_H
#define FIXED_POINT_PROFILE_H

#include "BeamProfile.h"
#include "TString.h"

class FixedPointProfile : public BeamProfile {

public:
	FixedPointProfile()
	{
		parX = 0.;
		parY = 0.;
	}

	std::pair<double, double> Sample() override {
		return {0., 0.};
	}

	TString ToString() const override{
		return "FixedPointProfile";
	}
};

#endif //FIXED_POINT_PROFILE_H
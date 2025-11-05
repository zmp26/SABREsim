#ifndef GAUSSIANPROFILEX_FIXEDPOINTY_H
#define GAUSSIANPROFILEX_FIXEDPOINTY_H

#include "BeamProfile.h"

class GaussianProfileX_FixedPointY : public BeamProfile{
private:
	std::mt19937 rng;
	std::normal_distribution<> distX;

public:
	GaussianProfileX_FixedPointY(double sigmaX, unsigned int seed = std::random_device{}())
		: rng(seed),
		distX(0.,sigmaX)
	{
		parX = sigmaX;
		parY = 0.;
	}

	std::pair<double,double> Sample() override {
		return {distX(rng), 0.};
	}

	TString ToString() const override {
		return "GaussianProfileX_FixedPointY";
	}
};

#endif//GAUSSIANPROFILEX_FIXEDPOINTY_H
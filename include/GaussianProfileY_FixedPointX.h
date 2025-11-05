#ifndef GAUSSIANPROFILEY_FIXEDPOINTX_H
#define GAUSSIANPROFILEY_FIXEDPOINTX_H

#include "BeamProfile.h"

class GaussianProfileY_FixedPointX : public BeamProfile{
private:
	std::mt19937 rng;
	std::normal_distribution<> distY;

public:
	GaussianProfileY_FixedPointX(double sigmaY, unsigned int seed = std::random_device{}())
		: rng(seed),
		distY(0.,sigmaY)
	{
		parX = 0.;
		parY = sigmaY;
	}

	std::pair<double,double> Sample() override {
		return {0.,distY(rng)};
	}

	TString ToString() const override {
		return "GaussianProfileY_FixedPointX";
	}
};

#endif//GAUSSIANPROFILEY_FIXEDPOINTX_H
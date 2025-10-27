#ifndef UNIFORM_PROFILE_H
#define UNIFORM_PROFILE_H

#include "BeamProfile.h"

class UniformProfile : public BeamProfile {
private:
	std::mt19937 rng;
	std::uniform_real_distribution<> distX;
	std::uniform_real_distribution<> distY;

public:
	UniformProfile(double xmax, double ymax, unsigned int seed = std::random_device{}())
		: rng(seed),
		  distX(-xmax, xmax),
		  distY(-ymax,ymax)
	{
		parX = xmax;
		parY = ymax;
	}

	std::pair<double, double> Sample() override {
		return {distX(rng), distY(rng)};
	}

	TString ToString() const override{
		return "UniformProfile";
	}
};

#endif //UNIFORM_PROFILE_H
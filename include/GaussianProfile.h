#ifndef GAUSSIAN_PROFILE_H
#define GAUSSIAN_PROFILE_H

#include "BeamProfile.h"

class GaussianProfile : public BeamProfile {
private:
	std::mt19937 rng;
	std::normal_distribution<> distX;
	std::normal_distribution<> distY;

public:
	GaussianProfile(double sigmaX, double sigmaY, unsigned int seed = std::random_device{}())
		: rng(seed),
		  distX(0.,sigmaX),
		  distY(0.,sigmaY)
	{
		parX = sigmaX;
		parY = sigmaY;
	}

	std::pair<double, double> Sample() override {
		return {distX(rng), distY(rng)};
	}

	TString ToString() const override{
		return "GaussianProfile";
	}

};

#endif //GAUSSIAN_PROFILE_H
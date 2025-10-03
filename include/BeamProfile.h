#ifndef BEAM_PROFILE_H
#define BEAM_PROFILE_H

#include <utility>
#include <random>

class BeamProfile{
public:
	virtual ~BeamProfile() = default;

	virtual std::pair<double, double> Sample() = 0;
};

#endif //BEAM_PROFILE_H
#ifndef TARGETANGULARSTRAGGLER
#define TARGETANGULARSTRAGGLER

#include "Distribution.h"
#include <memory>
#include <vector>
#include <random>

class TargetAngularStraggler {
public:
	explicit TargetAngularStraggler(std::unique_ptr<Distribution> dist);

	//sample a single value using distribution's internal rng
	double Sample();
	//sample a single value using some provided external rng
	double Sample(std::mt19937_64& rng);

	//sample single value of phi using this class's internal rng
	double SamplePhi();
	//sample single value of phi using some provided external rng
	double SamplePhi(std::mt19937_64& rng);

	//sample multiple values:
	std::vector<double> Sample(size_t n);
	std::vector<double> Sample(size_t n, std::mt19937_64& rng);

	//evaluate PDF at a point
	double pdf(double x) const;

private:
	std::unique_ptr<Distribution> distribution;

	std::mt19937_64 phi_rng_;
	std::uniform_real_distribution<double> phi_dist_{0., 360.};

};

#endif//TARGETANGULARSTRAGGLER
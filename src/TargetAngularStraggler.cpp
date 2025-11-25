#include "TargetAngularStraggler.h"
#include <random>

TargetAngularStraggler::TargetAngularStraggler(std::unique_ptr<Distribution> dist)
	: distribution(std::move(dist)), phi_rng_(std::random_device{}()) {}

//sample single value
double TargetAngularStraggler::Sample(){
	return distribution->Sample();
}

double TargetAngularStraggler::Sample(std::mt19937_64& rng){
	return distribution->Sample(rng);
}

//sample phi
double TargetAngularStraggler::SamplePhi(){
	return phi_dist_(phi_rng_);
}

double TargetAngularStraggler::SamplePhi(std::mt19937_64& rng){
	return phi_dist_(rng);
}

//sample multiple values:
std::vector<double> TargetAngularStraggler::Sample(size_t n){
	std::vector<double> samples(n);
	for(size_t i=0; i<n; i++){
		samples[i] = distribution->Sample();
	}
	return samples;
}

std::vector<double> TargetAngularStraggler::Sample(size_t n, std::mt19937_64& rng){
	std::vector<double> samples(n);
	for(size_t i=0; i<n; i++){
		samples[i] = distribution->Sample(rng);
	}
	return samples;
}

//evaluate PDF at point:
double TargetAngularStraggler::pdf(double x) const {
	return distribution->pdf(x);
}
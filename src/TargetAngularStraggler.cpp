#include "TargetAngularStraggler.h"
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

static inline std::string trim(const std::string& s){
	size_t start = s.find_first_not_of(" \t\n\r");
	size_t end = s.find_last_not_of(" \t\n\r");
	return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

TargetAngularStraggler::TargetAngularStraggler(std::unique_ptr<Distribution> dist)
	: distribution(std::move(dist)), mu_(0), sigma_(0), lambda_(0), phi_rng_(std::random_device{}()) {}

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

TargetAngularStraggler* TargetAngularStraggler::LoadFromConfigFile(const std::string& filename){
	std::ifstream file(filename);
	if(!file.is_open()){
		std::cerr << "Error: cannot open config file " << filename << std::endl;
		return nullptr; 
	}

	std::string line;
	double mu = 0., sigma = 0., lambda = 0.;
	while(std::getline(file,line)){
		line = trim(line);
		if(line.empty() || line[0]=='#') continue;

		auto pos = line.find('=');
		if(pos == std::string::npos) continue;

		std::string key = trim(line.substr(0,pos));
		std::string value = trim(line.substr(pos+1));

		if(key == "mu"){
			std::istringstream iss(value);
			iss >> mu;
		} else if(key == "sigma"){
			std::istringstream iss(value);
			iss >> sigma;
		} else if(key == "lambda"){
			std::istringstream iss(value);
			iss >> lambda;
		} else {
			std::cerr << "Warning: unknown key '" << key << "' in config file " << filename << "\n";
		}
	}

	TargetAngularStraggler* tas = new TargetAngularStraggler(std::make_unique<ExponentiallyModifiedGaussian>(mu,sigma,lambda));

	tas->SetMu_(mu);
	tas->SetSigma_(sigma);
	tas->SetLambda_(lambda);

	return tas;
}

double TargetAngularStraggler::GetDistMu(){
	return distribution->GetMu();
}

double TargetAngularStraggler::GetDistSigma(){
	return distribution->GetSigma();
}

double TargetAngularStraggler::GetDistLambda(){
	return distribution->GetLambda();
}

void TargetAngularStraggler::SetMu_(double mu){
	mu_ = mu;
}

void TargetAngularStraggler::SetSigma_(double sigma){
	sigma_ = sigma;
}

void TargetAngularStraggler::SetLambda_(double lambda){
	lambda_ = lambda;
}
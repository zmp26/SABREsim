#include <random>
#include <cmath>

class Distribution{
public:
	virtual ~Distribution() = default;

	//get single sample from distribution
	virtual double Sample(std::mt19937_64& extrng) = 0;//provide external rng
	virtual double Sample() = 0;//use internal rng

	//evaluate PDF for diagnostics
	virtual double pdf(double x) const = 0;
};



class ExponentiallyModifiedGaussian : public Distribution {
public:
	ExponentiallyModifiedGaussian(double mu, double sigma, double lambda, uint64_t seed = std::random_device{}())
		: mu(mu), sigma(sigma), lambda(lambda), normal(mu, sigma), exponential(lambda), rng(seed) {}

	double Sample(std::mt19937_64& extrng) override {
		return normal(extrng) + exponential(extrng);
	}

	double Sample() override {
		return normal(rng) + exponential(rng);
	}

	double pdf(double x) const override {

		double ls2 = lambda * sigma * sigma;

		double z = (lambda * (mu + ls2 - x));
		double arg = (x - mu - ls2) / (std::sqrt(2.0)*sigma);
		return (lambda / 2.0) * std::exp(z + ls2/2.0) * std::erfc(arg);
	}

private:
	double mu, sigma, lambda;
	std::normal_distribution<double> normal;
	std::exponential_distribution<double> exponential;
	std::mt19937_64 rng;
};
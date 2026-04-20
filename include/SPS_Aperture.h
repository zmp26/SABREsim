#ifndef SPS_APERTURE_H
#define SPS_APERTURE_H

#include <iostream>
#include <cmath>
#include <random>
#include "Vec3.h"

class SPS_Aperture{
public:
	SPS_Aperture(double r, 
				 double thetaMin, double thetaMax,
				 double phiMin, double phiMax,
				 double sigE=0., double sigTheta=0., double sigPhi=0.);

	~SPS_Aperture();

	bool IsDetected(const Vec3& traj, const Vec3& origin) const;

	double GetSmearedEnergy(double energy) const;
	double GetSmearedTheta(double theta) const;
	double GetSmearedPhi(double phi) const;

	void SetEnergyResolution(double sigma) { sigmaE = sigma; }
	void SetThetaResolution(double sigma) { sigmaTheta = sigma; }
	void SetPhiResolution(double sigma) { sigmaPhi = sigma; }

private:
	const double DEGRAD = M_PI / 180.;
	const double RADDEG = 180. / M_PI;

	double ApplyGaussian(double value, double sigma) const;


	double radius, thetaMin, thetaMax, phiMin, phiMax;

	double sigmaE, sigmaTheta, sigmaPhi;

	mutable std::mt19937 gen;
};

#endif//SPS_APERTURE_H
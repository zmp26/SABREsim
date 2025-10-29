#ifndef TARGET_ENERGY_LOSS_H
#define TARGET_ENERGY_LOSS_H

#include "TF1.h"
#include <string>
#include <vector>

class TargetEnergyLoss{
public:
	TargetEnergyLoss(const std::string& funcStr,
					 const std::vector<double>& params,
					 double arealDensity,					//in ug/cm^2
					 double materialDensity,				//in g/cm^3
					 const std::string& tostringmsg);

	~TargetEnergyLoss();

	TargetEnergyLoss(const TargetEnergyLoss&) = delete;
	TargetEnergyLoss& operator=(const TargetEnergyLoss&) = delete;

	double ApplyEnergyLoss(double energy_MeV, double theta_deg);

	double GetPathLength(double theta_deg) const;

	void SetLossFunction(const std::string& funcStr, const std::vector<double>& params);

	double EvaluateLossFunction(double energy_MeV) const;

	static TargetEnergyLoss* LoadFromConfigFile(const std::string& filename);

	std::string ToString() const { return ToStringMessage; }

private:
	TF1* lossFunction;										//ROOT function to describe energy loss (keV) per unit areal density at initial energy E
	double arealDensity;									//in ug/cm^2
	double materialDensity;									//in g/cm^3
	double linearThickness;									//cm, computed from arealDensity / material Density
	std::string ToStringMessage;
};

#endif

/*

//old version (constant subtraction)
#ifndef TARGET_ENERGY_LOSS_H
#define TARGET_ENERGY_LOSS_H

class TargetEnergyLoss{
public:
	TargetEnergyLoss();

	double ApplyEnergyLoss(double energyMeV) const;

private:
	static constexpr double kEnergyLossKeV = 280.0; //fixed energy loss for 6Li ions in 75 ug/cm2 at 45 degrees - will need to change this later!!!!!!!!

};

#endif

*/
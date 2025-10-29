#ifndef SABRE_DEADLAYERMODEL_H
#define SABRE_DEADLAYERMODEL_H

#include "TF1.h"
#include "Vec3.h"
#include <string>
#include <vector>

class SABRE_DeadLayerModel{
public:
	SABRE_DeadLayerModel(const std::string& funcStr,
						 const std::vector<double>& params,
						 const std::string& tostringmsg);

	~SABRE_DeadLayerModel();

	SABRE_DeadLayerModel(const SABRE_DeadLayerModel&) = delete;
	SABRE_DeadLayerModel& operator=(const SABRE_DeadLayerModel&) = delete;

	double ApplyEnergyLoss(double energy_MeV, Vec3& trajectory, Vec3& detectorNormal);

	double GetPathLength(double theta_deg) const;

	void SetLossFunction(const std::string& funcStr, const std::vector<double>& params);

	double EvaluateLossFunction(double energy_MeV) const;

	static SABRE_DeadLayerModel* LoadFromConfigFile(const std::string& filename);

	std::string ToString() const { return ToStringMessage; }

private:
	TF1* lossFunction;
	const double materialDensity = 2.33;		// g/cm^3
	const double linearThickness = 5e-6;		// cm			(50nm deadlayer = 5e-6 cm)
	const double arealDensity = materialDensity*linearThickness;					// g/cm^2
	std::string ToStringMessage;
};


#endif



//previous:
/*
#ifndef SABRE_DEADLAYERMODEL_H
#define SABRE_DEADLAYERMODEL_H

class SABRE_DeadLayerModel{
public:
	SABRE_DeadLayerModel();

	double ApplyEnergyLoss(double energyMeV) const;

private:
	static constexpr double kDeadLayerEnergyLossKeV = 25.0;
};

#endif
*/
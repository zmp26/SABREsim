#include "TargetEnergyLoss.h"

TargetEnergyLoss::TargetEnergyLoss() {}

double TargetEnergyLoss::ApplyEnergyLoss(double energyMeV) const {
	double energyLossMeV = kEnergyLossKeV / 1000.;
	double energyOut = energyMeV - energyLossMeV;

	if(energyOut < 0){
		energyOut = 0.;
	}

	return energyOut;
}
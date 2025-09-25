#include "SABRE_DeadLayerModel.h"

SABRE_DeadLayerModel::SABRE_DeadLayerModel() {}

double SABRE_DeadLayerModel::ApplyEnergyLoss(double energyMeV) const {
	double energyLossMeV = kDeadLayerEnergyLossKeV / 1000.;
	double energyOut = energyMeV - energyLossMeV;

	if(energyOut < 0.){
		energyOut = 0.;
	}

	return energyOut;
}
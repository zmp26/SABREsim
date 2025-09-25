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
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
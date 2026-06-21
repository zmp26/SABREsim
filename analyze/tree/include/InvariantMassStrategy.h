#ifndef INVARIANTMASSSTRATEGY_H
#define INVARIANTMASSSTRATEGY_H

#include "ReconstructionStrategy.h"

class InvariantMassStrategy : public ReconstructionStrategy {
public:
	void ResolveMomentum(KinematicNode* root) const override;
	double CalculateExcitation(const KinematicNode* node) const override;
	double CalculateDecayEnergy(const KinematicNode* node) const override;
};

#endif//INVARIANTMASSSTRATEGY_H
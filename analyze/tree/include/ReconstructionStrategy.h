#ifndef RECONSTRUCTIONSTRATEGY_H
#define RECONSTRUCTIONSTRATEGY_H

#include "KinematicNode.h"

class ReconstructionStrategy {
public:
	virtual ~ReconstructionStrategy() = default;

	virtual void ResolveMomentum(KinematicNode* root) const = 0;

	virtual double CalculateExcitation(const KinematicNode* node) const = 0;

	virtual double CalculateDecayEnergy(const KinematicNode* node) const = 0;
};

#endif//RECONSTRUCTIONSTRATEGY_H
#ifndef MISSINGMASSSTRATEGY_H
#define MISSINGMASSSTRATEGY_H

#include "ReconstructionStrategy.h"
#include "TLorentzVector.h"

class MissingMassStrategy : public ReconstructionStrategy {
private:
	TLorentzVector beam, target;
	std::string missingNodeName;

	void AssignMissingMomentum(KinematicNode* node, const TLorentzVector& p) const;

public:
	MissingMassStrategy(TLorentzVector b, TLorentzVector t) : beam(b), target(t) {};

	void ResolveMomentum(KinematicNode* root) const override;
	double CalculateExcitation(const KinematicNode* node) const override;
	double CalculateDecayEnergy(const KinematicNode* node) const override;
};

#endif//MISSINGMASSSTRATEGY_H
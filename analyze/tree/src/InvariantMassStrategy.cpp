#include "InvariantMassStrategy.h"

void InvariantMassStrategy::ResolveMomentum(KinematicNode* root) const {
	root->CalculateMomentum();
}

double InvariantMassStrategy::CalculateExcitation(const KinematicNode* node) const {
	return node->p.M() - node->restMass;
}

double InvariantMassStrategy::CalculateDecayEnergy(const KinematicNode* node) const {
	if(node->children.empty()) return 0.;

	double sumMass = 0.;
	for(const auto& child : node->children){
		sumMass += child->restMass;
	}

	return node->p.M() - sumMass;
}
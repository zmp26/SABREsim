#include "MissingMassStrategy.h"

void MissingMassStrategy::ResolveMomentum(KinematicNode* root) const {

	TLorentzVector totalSystem = beam + target;

	TLorentzVector detectedSum;

	TLorentzVector missingP = totalSystem - detectedSum;


	root->CalculateMomentum();

}
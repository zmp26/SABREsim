#include "EventMapper.h"
#include <stdexcept>

bool EventMapper::Initialize(const std::string& config){
	return tree.LoadTree(config);
}

void EventMapper::MapData(const std::vector<int>& perm, const std::vector<TLorentzVector>& hits){
	auto leaves = tree.GetLeafNodes();

	if(perm.size() != leaves.size()){
		throw std::runtime_error("Permutation size mismatch!");
	}

	for(size_t i=0; i<leaves.size(); i++){
		leaves[i]->p = hits[perm[i]];
		leaves[i]->is_leaef = true;
	}
}

void EventMapper::ResetTree(){
	tree.ResetAllNodes();
}

void EventMapper::ProcessEvent(){
	if(tree.root){
		tree.root->CalculateMomentum();
	}
}

KinematicNode* EventMapper::GetRoot(){
	return tree.root.get();
}

//void EventMapper::ResolveMissingMass(const TLorentzVector& beam, const TLorentzVector& target, const TLorentzVector& ejectile){}?
#ifndef EVENTMAPPER_H
#define EVENTMAPPER_H

#include "DecayTree.h"
#include "ReconstructionEngine.h"
#include <vector>
#include <string>

class EventMapper {
private:
	DecayTree tree;
	ReconstructionEngine engine;

public:

	EventMapper() = default;

	//loads specific decay tree topology from config file
	bool Initialize(const std::string& config);//{ return tree.LoadTree(config);}

	// void MapData(const std::string& nodeName, const TLorentzVector& p){
	// 	KinematicNode* node = tree.GetNode(nodeName);
	// 	if(node){
	// 		node->p = p;
	// 		node->is_leaf = true;
	// 	}
	// }

	//maps TLorentzVectors to DecayTree leaves/KinematicNodes for a given permutation perm
	void MapData(const std::vector<int>& perm, const std::vector<TLorentzVector>& hits);//{

	// 	std::vector<std::string> leafNames = {"Frag1", "Frag2", "Frag3"};

	// 	for(size_t i=0; i<leafNames.size(); i++){
	// 		KinematicNode* node = tree.FindNode(leafNames[i]);
	// 		if(node){
	// 			node->p=hits[perm[i]];
	// 			node->is_leaf = true;
	// 		}
	// 	}

	// }

	//Resets the momentum in all nodes of the tree to 0
	void ResetTree();

	//Triggers calculation of event
	void ProcessEvent();//{
	// 	if(tree.root){
	// 		tree.root->CalculateMomentum();
	// 	}
	// }

	//Getter for DecayTree root node (reaction recoil)
	KinematicNode* GetRoot();// { return tree.root.get(); }

	// void ResolveMissingMass(const TLorentzVector& beam, const TLorentzVector& target, const TLorentzVector& ejectile){
	// 	TLorentzVector pRecoil = (beam + target) - ejectile;

	// 	KinematicNode* recoilNode = tree.FindNode("Recoil");
	// 	if(recoilNode){
	// 		recoilNode->p = pRecoil;
	// 		recoilNode->is_leaf = true;
	// 		recoilNode->is_missing = true;
	// 	}

	// 	ProcessEvent();
	// }

	//handles special case of missing particles
	void ResolveMissingMass(const TLorentzVector& beam, const TLorentzVector& target, const TLorentzVector& ejectile);

};


#endif//EVENTMAPPER_H
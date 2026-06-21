#include "DecayTree.h"
#include <fstream>
#include <sstream>
#include <iostream>

bool DecayTree::LoadTree(const std::string& filename){
	std::ifstream file(filename);
	std::string line;

	while(std::getline(file, line)){
		std::stringstream ss(line);
		std::string name, massStr, parentName;

		std::getline(ss, name, ',');
		std::getline(ss, massStr, ',');
		std::getline(ss, parentName, ',');

		double mass = std::stod(massStr);
		auto newNode = std::make_unique<KinematicNode>(name, mass);

		KinematicNode* ptr = newNode.get();
		nodeRegistry[name] = ptr;

		if(parentName == "ROOT"){
			root = std::move(newNode);
		} else {
			//nodeRegistry[parentName]->AddChild(std::move(newNode));
			auto it = nodeRegistry.find(parentName);
			if(it != nodeRegistry.end()){
				it->second->AddChild(std::move(newNode));
			} else {
				std::cerr << "Error: Parent '" << parentName << "' not found for node '" << name << "'" << std::endl;
				return false;
			}
		}
	}
	return true;
}

std::unique_ptr<DecayTree> DecayTree::CreateSequential4Body(double mRecoil, double mInt, double mFrag1, double mFrag2, double mFrag3){
	auto tree = std::make_unique<DecayTree>();

	auto rootNode = std::make_unique<KinematicNode>("Recoil",mRecoil);
	auto frag1 = std::make_unique<KinematicNode>("Frag1",mFrag1);
	auto intermediate = std::make_unique<KinematicNode>("Intermediate",mInt);
	auto frag2 = std::make_unique<KinematicNode>("Frag2",mFrag2);
	auto frag3 = std::make_unique<KinematicNode>("Frag3",mFrag3);

	//set leaf status
	frag1->is_leaf = true;
	frag2->is_leaf = true;
	frag3->is_leaf = true;


	//link nodes together to form tree
	intermediate->AddChild(std::move(frag2));
	intermediate->AddChild(std::move(frag3));
	rootNode->AddChild(std::move(frag1));
	rootNode->AddChild(std::move(intermediate));

	tree->root = std::move(rootNode);
	return tree;
}

std::vector<KinematicNode*> DecayTree::GetLeafNodes(){
	std::vector<KinematicNode*> leaves;
	for(auto const& [name, node] : nodeRegistry){
		if(node->children.empty()) leaves.push_back(node);
	}
	return leaves;
}

void DecayTree::ResetAllNodes(){
	for(auto const& [name, node] : nodeRegistry){
		node->Reset();
	}
}

KinematicNode* DecayTree::FindNode(const std::string& name){
	auto it = nodeRegistry.find(name);
	if(it != nodeRegistry.end()){
		return it->second;
	}
	return nullptr;
}
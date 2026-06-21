#ifndef KINEMATICNODE_H
#define KINEMATICNODE_H

#include <string>
#include <vector>
#include <memory>
#include "TLorentzVector.h"

struct KinematicNode {

	std::string name;
	double restMass;
	TLorentzVector p;
	bool is_leaf = false;
	bool is_missing = false;

	std::vector<std::unique_ptr<KinematicNode>> children;

	KinematicNode(std::string n, double m) : name(n), restMass(m) {}

	void AddChild(std::unique_ptr<KinematicNode> child){
		children.push_back(std::move(child));
	}

	//recursively calculate 4vector:
	//if leaf, return own p, else return sum of children
	TLorentzVector CalculateMomentum(){
		if(is_leaf) return p;

		TLorentzVector sum;
		for(const auto& child : children) sum += child->CalculateMomentum();

		p = sum;
		return p;
	}

	void Reset(){
		p.SetPxPyPzE(0,0,0,0);
		is_leaf = false;
		for(auto& child : children) child->Reset();
	}

};

#endif//KINEMATICNODE_H
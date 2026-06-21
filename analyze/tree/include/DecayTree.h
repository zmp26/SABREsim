#ifndef DECAYTREE_H
#define DECAYTREE_H

#include "KinematicNode.h"
#include <memory>
#include <string>
#include <map>

class DecayTree {
private:
	std::map<std::string, KinematicNode*> nodeRegistry;
public:
	std::unique_ptr<KinematicNode> root;

	static std::unique_ptr<DecayTree> CreateSequential4Body(double mRecoil, double mInt, double mFrag1, double mFrag2, double mFrag3);

	bool LoadTree(const std::string& filename);

	KinematicNode* FindNode(const std::string& name);

	std::vector<KinematicNode*> GetLeafNodes();

	void ResetAllNodes();
};

#endif//DECAYTREE_H
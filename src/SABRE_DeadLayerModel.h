#ifndef SABRE_DEADLAYERMODEL_H
#define SABRE_DEADLAYERMODEL_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

#include "Vec3.h"

class SABRE_DeadLayerModel{
public:
	SABRE_DeadLayerModel();
	~SABRE_DeadLayerModel();
	void LoadSRIMData(const std::string& particleType, const std::string& filename);
	double ApplyDeadLayerCorrection(const std::string& particleType, double energy, double theta, double phi) const;

private:
	//store SRIM data as a lookup table
	std::map<std::string,std::vector<std::pair<double,double>>> fStoppingPower; //particle type -> (E,dE/dx)
	//std::map<std::string,std::vector<std::pair<double,std::array<double,3>>>> fStragglingData; //particle type -> (E,{range,longstrat,latstrat})
	double GetEnergyLoss(const std::string& particleType, double energy, double pathLength) const;
};

#endif
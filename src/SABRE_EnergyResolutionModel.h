#ifndef SABRE_ENERGYRESOLUTIONMODEL_H
#define SABRE_ENERGYRESOLUTIONMODEL_H

#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <string>

class SABRE_EnergyResolutionModel{
public:
	static constexpr int NUM_RINGS = 16;
	static constexpr int NUM_WEDGES = 8;


	SABRE_EnergyResolutionModel(double sigmaMeV=0.05, double thresholdMeV=1.0);

	//setters
	void setRingResolution(int ring, double sigmaMeV);
	void setWedgeResolution(int wedge, double sigmaMeV);
	void setRingThreshold(int ring, double thresholdMeV);
	void setWedgeThreshold(int wedge, double thresholdMeV);
	void setNormalDist(double sigma, double centroid=0.);

	//getters
	double getRingResolution(int ring) { if(ring>=0 && ring<NUM_RINGS) { return ringSigmas[ring]; } else { return -1.; } }
	double getWedgeResolution(int wedge){ if(wedge>=0 && wedge<NUM_WEDGES) { return wedgeSigmas[wedge]; } else { return -1; } }
	double getRingThreshold(int ring) { if(ring>=0 && ring<NUM_RINGS) { return ringThresholds[ring]; } else { return -1.; } }
	double getWedgeThreshold(int wedge) { if(wedge>=0 && wedge<NUM_WEDGES) { return wedgeThresholds[wedge]; } else { return -1.; } }
	std::normal_distribution<double> getRingNormalDist(int ring);
	std::normal_distribution<double> getWedgeNormalDist(int wedge);

	//load resolutions and thresholds from a file
	bool loadFromFile(const std::string& filename);

	//main interface
	bool detectEnergyInRing(int ring, double kinEnergyMeV, double& detectEnergyMeV);//detectEnergyMeV updated to kinEnergy with resolution applied if above threshold
	bool detectEnergyInWedge(int wedge, double kinEnergyMeV, double& detectEnergyMeV);//detectEnergyMeV updated to kinEnergy with resolution applied if above threshold

private:
	std::vector<double> ringSigmas;
	std::vector<double> wedgeSigmas;
	std::vector<double> ringThresholds;
	std::vector<double> wedgeThresholds;

	double applyEnergyResolution(double kinEnergyMeV, double sigma);
};

#endif
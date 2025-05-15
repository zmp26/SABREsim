#include "SABRE_EnergyResolutionModel.h"

SABRE_EnergyResolutionModel::SABRE_EnergyResolutionModel(double sigmaMeV, double thresholdMeV)
	: ringSigmas(NUM_RINGS, sigmaMeV),
	  wedgeSigmas(NUM_WEDGES, sigmaMeV),
	  ringThresholds(NUM_RINGS,thresholdMeV),
	  wedgeThresholds(NUM_WEDGES,thresholdMeV)
{
	//nothing
}

std::normal_distribution<double> SABRE_EnergyResolutionModel::getRingNormalDist(int ring){
		if(ring>=0 && ring < NUM_RINGS){
			return std::normal_distribution<double>(0.,getRingResolution(ring));
		} else {
			return std::normal_distribution<double>(0.,0.);
		}
	}
std::normal_distribution<double> SABRE_EnergyResolutionModel::getWedgeNormalDist(int wedge){
		if(wedge>=0 && wedge < NUM_WEDGES){
			return std::normal_distribution<double>(0.,getWedgeResolution(wedge));
		} else {
			return std::normal_distribution<double>(0.,0.);
		}
	}

void SABRE_EnergyResolutionModel::setRingResolution(int ring, double sigmaMeV){
	if(ring>=0 && ring < NUM_RINGS) ringSigmas[ring] = sigmaMeV;
}

void SABRE_EnergyResolutionModel::setWedgeResolution(int wedge, double sigmaMeV){
	if(wedge>=0 && wedge<NUM_WEDGES) wedgeSigmas[wedge] = sigmaMeV;
}

void SABRE_EnergyResolutionModel::setRingThreshold(int ring, double thresholdMeV){
	if(ring>=0 && ring<NUM_RINGS) ringThresholds[ring] = thresholdMeV;
}

void SABRE_EnergyResolutionModel::setWedgeThreshold(int wedge, double thresholdMeV){
	if(wedge>=0 && wedge<NUM_WEDGES) wedgeThresholds[wedge] = thresholdMeV;
}

double SABRE_EnergyResolutionModel::applyEnergyResolution(double kinEnergyMeV, double sigma){
	std::random_device rd;
	std::mt19937 gen(rd());//mersenne twister rng
	std::normal_distribution<double> dist(0.,sigma);
	double sample = dist(gen);

	return kinEnergyMeV + sample;
}

bool SABRE_EnergyResolutionModel::detectEnergyInRing(int ring, double kinEnergyMeV, double& detectedEnergyMeV){
	if(ring<=0 || ring >= NUM_RINGS) return false;

	detectedEnergyMeV = applyEnergyResolution(kinEnergyMeV,getRingResolution(ring));
	return detectedEnergyMeV >= getRingThreshold(ring);
}

bool SABRE_EnergyResolutionModel::detectEnergyInWedge(int wedge, double kinEnergyMeV, double& detectedEnergyMeV){
	if(wedge<=0 || wedge >= NUM_WEDGES) return false;

	detectedEnergyMeV = applyEnergyResolution(kinEnergyMeV,getWedgeResolution(wedge));
	return detectedEnergyMeV >= getWedgeThreshold(wedge);
}

bool SABRE_EnergyResolutionModel::loadFromFile(const std::string& filename){
	std::ifstream infile(filename);
	if(!infile.is_open()) return false;

	int ring, wedge;
	double sigma, threshold;

	while(infile >> ring >> wedge >> sigma >> threshold){
		if(ring>=0 && ring<NUM_RINGS){
			ringSigmas[ring] = sigma;
			ringThresholds[ring] = threshold;
		}
		if(wedge>=0 && wedge<NUM_WEDGES){
			wedgeSigmas[wedge] = sigma;
			wedgeThresholds[wedge] = threshold;
		}
	}

	infile.close();
	return true;
}
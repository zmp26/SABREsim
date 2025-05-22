#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

#include "SABRE_DeadLayerModel.h"

SABRE_DeadLayerModel::SABRE_DeadLayerModel(){

}

SABRE_DeadLayerModel::~SABRE_DeadLayerModel(){

}

void SABRE_DeadLayerModel::LoadSRIMData(const std::string& particleType, const std::string& filename){
	std::ifstream infile(filename);
	if(!infile.is_open()){
		throw std::runtime_error("Failed to open SRIM file " + filename); 
	}

	std::string line;
	bool dataSectionStarted = false;
	std::vector<std::pair<double,double>> data;
	//std::vector<std::pair<double,double[3]>> stragglingdata;

	while(std::getline(infile,line)){
		if(!dataSectionStarted){
			if(line.find("        Ion        dE/dx") != std::string::npos){
				dataSectionStarted = true;
				std::getline(infile,line);//skip next line
			}
			continue;
		}

		if(line.empty() || line.find("---") != std::string::npos) break;

		std::istringstream iss(line);
		std::string energyStr, unit;
		double energy, dEdxElec, dEdxNucl, range, longstrag, latstrag;

		iss >> energyStr;

		if(energyStr.find("MeV") != std::string::npos){
			energyStr.erase(energyStr.find("MeV"),3);
			energy = std::stod(energyStr);
		} else if(energyStr.find("keV") != std::string::npos){
			energyStr.erase(energyStr.find("keV"),3);
			energy = std::stod(energyStr)/1000.;
		} else {
			energy = std::stod(energyStr)/1000.;
		}

		iss >> dEdxElec >> dEdxNucl >> range >> longstrag >> latstrag;
		double totalStopping = dEdxElec + dEdxNucl;

		data.emplace_back(energy,totalStopping);
		//stragglingdata.emplace_back(energy,{range,longstrag,latstrag});
	}

	if(data.empty()){
		throw std::runtime_error("No data found in SRIM file: " + filename);
	}

	fStoppingPower[particleType] = std::move(data);
}
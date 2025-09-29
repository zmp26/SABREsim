#include "SABRE_DeadLayerModel.h"
#include "Vec3.h"
#include "ConsoleColorizer.h"
#include "TMath.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

SABRE_DeadLayerModel::SABRE_DeadLayerModel(const std::string& funcStr,
										   const std::vector<double>& params) {
	lossFunction = new TF1("DeadLayerEnergyLossFunction", funcStr.c_str(), 0, 10000);
	for(size_t i=0; i<params.size(); i++){
		lossFunction->SetParameter(i,params[i]);
	}
}

SABRE_DeadLayerModel::~SABRE_DeadLayerModel(){
	delete lossFunction;
}

static inline std::string trim(const std::string& s){
	size_t start = s.find_first_not_of(" \t\n\r");
	size_t end = s.find_last_not_of(" \t\n\r");
	return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

SABRE_DeadLayerModel* SABRE_DeadLayerModel::LoadFromConfigFile(const std::string& filename){
	std::ifstream file(filename);
	if(!file.is_open()){
		std::cerr << "Error: cannot open config file " << filename << std::endl;
		return nullptr;
	}

	std::string funcStr;
	std::vector<double> params;

	std::string line;
	while(std::getline(file,line)){
		line = trim(line);
		if(line.empty() || line[0]=='#') continue;

		auto pos = line.find('=');
		if(pos == std::string::npos) continue;

		std::string key = trim(line.substr(0,pos));
		std::string value = trim(line.substr(pos+1));

		if(key=="funcStr"){
			if(!value.empty() && value.front() == '"' && value.back() == '"'){
				value = value.substr(1,value.size()-2);
			}
			funcStr = value;
		} else if(key=="params"){
			params.clear();
			std::istringstream iss(value);
			double val;
			while(iss >> val){
				params.push_back(val);
			}
		} else {
			std::cerr << "Warning: unknown key '" << key << "' in config file\n";
		}
	}

	if(funcStr.empty()){
		ConsoleColorizer::PrintRed("Error: funcStr not specified in config file\n");
		return nullptr;
	}

	if(params.empty()){
		ConsoleColorizer::PrintRed("Error: funcStr not specified in config file\n");
		return nullptr;
	}

	return new SABRE_DeadLayerModel(funcStr,params);
}

double SABRE_DeadLayerModel::GetPathLength(double theta_deg) const {
	double theta_rad = TMath::DegToRad()*theta_deg;
	return fabs(linearThickness/std::cos(theta_rad));
}

double SABRE_DeadLayerModel::EvaluateLossFunction(double energy_MeV) const {
	if(!lossFunction) throw std::runtime_error("Loss function not initialized");
	return lossFunction->Eval(energy_MeV*1000.)/1000.;//MeV/ug/cm^2
}

double SABRE_DeadLayerModel::ApplyEnergyLoss(double energy_MeV, Vec3& trajectory, Vec3& detectorNormal){
	//get angle between detector normal and particle trajectory:
	double angle_deg = TMath::RadToDeg()*(TMath::ACos(detectorNormal.Dot(trajectory)));///(trajectory.Mag()*detectorNormal.Mag())));
	double path_cm = GetPathLength(angle_deg);

	double effectiveArealDensity = path_cm*materialDensity*1e6; //ug/cm^2

	double dEdx_MeV = EvaluateLossFunction(energy_MeV); //MeV/ug/cm^2

	double deltaE_MeV = (dEdx_MeV*effectiveArealDensity); //MeV

	double energy_out_MeV = energy_MeV - (deltaE_MeV);

	//std::cout << "path_cm = " << path_cm << std::endl;

	if(energy_out_MeV < 0) energy_out_MeV = 0.;

	return energy_out_MeV;
}

void SABRE_DeadLayerModel::SetLossFunction(const std::string& funcStr, const std::vector<double>& params){
	delete lossFunction;
	lossFunction = new TF1("DeadLayerEnergyLossFunction", funcStr.c_str(),0,20000);
	for(size_t i=0; i<params.size(); i++){
		lossFunction->SetParameter(i,params[i]);
	}
}

// previous:
/*
#include "SABRE_DeadLayerModel.h"

SABRE_DeadLayerModel::SABRE_DeadLayerModel() {}

double SABRE_DeadLayerModel::ApplyEnergyLoss(double energyMeV) const {
	double energyLossMeV = kDeadLayerEnergyLossKeV / 1000.;
	double energyOut = energyMeV - energyLossMeV;

	if(energyOut < 0.){
		energyOut = 0.;
	}

	return energyOut;
}
*/
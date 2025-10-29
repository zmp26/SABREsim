#include "TargetEnergyLoss.h"
#include "TMath.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

TargetEnergyLoss::TargetEnergyLoss(const std::string& funcStr,
								   const std::vector<double>& params,
								   double arealDensity_,
								   double materialDensity_,
								   const std::string& tostringmsg)
	: arealDensity(arealDensity_), materialDensity(materialDensity_), ToStringMessage(tostringmsg) {
	if(materialDensity <= 0) throw std::invalid_argument("Material density must be > 0");
	if(arealDensity <=0) throw std::invalid_argument("Areal density must be > 0");

	double arealdensity_gcm2 = arealDensity*1e-6;
	linearThickness = arealdensity_gcm2 / materialDensity;//in cm

	lossFunction = new TF1("TargetEnergyLossFunction", funcStr.c_str(), 0, 20000);
	for(size_t i=0; i<params.size(); i++){
		lossFunction->SetParameter(i,params[i]);
	}
}

TargetEnergyLoss::~TargetEnergyLoss(){
	delete lossFunction;
}

static inline std::string trim(const std::string& s){
	size_t start = s.find_first_not_of(" \t\n\r");
	size_t end = s.find_last_not_of(" \t\n\r");
	return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

TargetEnergyLoss* TargetEnergyLoss::LoadFromConfigFile(const std::string& filename){
	std::ifstream file(filename);
	if(!file.is_open()){
		std::cerr << "Error: Cannot open config file " << filename << std::endl;
		return nullptr;
	}

	std::string funcStr;
	std::vector<double> params;
	double arealDensity_=-1;
	double materialDensity_=-1;

	std::string line;
	std::string tostringmsg;
	while(std::getline(file,line)){
		line = trim(line);
		if(line.empty() || line[0]=='#') continue;

		auto pos = line.find('=');
		if(pos == std::string::npos) continue;

		std::string key = trim(line.substr(0,pos));
		std::string value = trim(line.substr(pos+1));

		if(key == "funcStr"){
			if(!value.empty() && value.front() == '"' && value.back() == '"' ){
				value = value.substr(1,value.size()-2);
			}
			funcStr = value;
		} else if(key == "params"){
			params.clear();
			std::istringstream iss(value);
			double val;
			while(iss >> val){
				params.push_back(val);
			}
		} else if(key == "arealDensity"){
			try{
				arealDensity_ = std::stod(value);
			} catch(...){
				std::cerr << "Error parsing arealDensity\n";
				return nullptr;
			}
		} else if(key == "materialDensity"){
			try{
				materialDensity_ = std::stod(value);
			} catch (...){
				std::cerr << "Error parsing materialDensity\n";
				return nullptr;
			}
		} else if(key == "ToString"){
			if(!value.empty() && value.front() == '"' && value.back() == '"'){
				value = value.substr(1,value.size()-2);
			}
			tostringmsg = value;
		} else {
			std::cerr << "Warning: unknown key '" << key << "' in config file\n";
		}
	}

	if(funcStr.empty()){
		std::cerr << "Error: funcStr not specified in config file\n";
		return nullptr;
	}

	if(params.empty()){
		std::cerr << "Error: params not specified in config file\n";
		return nullptr;
	}

	if(arealDensity_ <= 0){
		std::cerr << "Error: Invalid areal density\n";
		return nullptr;
	}

	if(materialDensity_ <= 0){
		std::cerr << "Error: Invalid material density\n";
		return nullptr;
	}

	return new TargetEnergyLoss(funcStr, params, arealDensity_, materialDensity_, tostringmsg);
}

double TargetEnergyLoss::GetPathLength(double theta_deg) const {
	double theta_rad = TMath::DegToRad()*theta_deg;
	return fabs(linearThickness/std::cos(theta_rad));
}

double TargetEnergyLoss::EvaluateLossFunction(double energy_MeV) const {
	if(!lossFunction) throw std::runtime_error("Loss function not initialized!");
	return lossFunction->Eval(energy_MeV*1000.)/1000.;//MeV/ug/cm^2
}

double TargetEnergyLoss::ApplyEnergyLoss(double energy_MeV, double theta_deg){
	double path_cm = GetPathLength(theta_deg);//cm

	//effective areal density along path length in ug/cm^2:
	double effectiveArealDensity = path_cm*materialDensity*1e6;// ug/cm^2

	double dEdx_MeV = EvaluateLossFunction(energy_MeV);// MeV/ug/cm^2

	double deltaE_MeV = (dEdx_MeV*effectiveArealDensity);// MeV

	double energy_out_MeV = energy_MeV - (deltaE_MeV);

	if(energy_out_MeV < 0) energy_out_MeV = 0.;

	return energy_out_MeV;
}

void TargetEnergyLoss::SetLossFunction(const std::string& funcStr, const std::vector<double>& params){
	delete lossFunction;
	lossFunction = new TF1("TargetEnergyLossFunction",funcStr.c_str(),0,20000);
	for(size_t i=0; i<params.size(); i++){
		lossFunction->SetParameter(i, params[i]);
	}
}

/*

//old version (constant subtraction)

#include "TargetEnergyLoss.h"

TargetEnergyLoss::TargetEnergyLoss() {}

double TargetEnergyLoss::ApplyEnergyLoss(double energyMeV) const {
	double energyLossMeV = kEnergyLossKeV / 1000.;
	double energyOut = energyMeV - energyLossMeV;

	if(energyOut < 0){
		energyOut = 0.;
	}

	return energyOut;
}


*/
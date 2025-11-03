#include "SimConfig.h"
#include "ConsoleColorizer.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "TString.h"

SimConfig::SimConfig(const std::string& filename)
	: filename_(filename),
	  detmc_version_(0),
	  beam_parX_(0.), beam_parY_(0.),
	  beam_offsetX_(0.), beam_offsetY_(0.),
	  beam_energy_(0.), recoil_excitation_energy_(0.)
{

	targetLoss_par_.resize(4, "none");
	deadLayerLoss_par_.resize(4, "none");

}

std::string SimConfig::Trim(const std::string& s){
	auto start = s.find_first_not_of(" \t\r\n");
	auto end = s.find_last_not_of(" \t\r\n");
	if(start == std::string::npos) return "";
	return s.substr(start, end-start+1);
}

bool SimConfig::Parse(){
	std::ifstream infile(filename_);
	if(!infile.is_open()){
		ConsoleColorizer::PrintRed(Form("Error! Could not open config file: %s",filename_.data()));
		return false;
	}

	std::string line, section;
	while(std::getline(infile, line)){
		line = Trim(line);
		if(line.empty() || line[0] == '#') continue;

		if(line.front() == '[' && line.back() == ']'){
			section = line.substr(1, line.size() -2);
			continue;
		}

		auto eq = line.find('=');
		if(eq == std::string::npos) continue;
		std::string key = Trim(line.substr(0,eq));
		std::string val = Trim(line.substr(eq+1));

		if(section == "General"){
			if(key == "detmc_version") detmc_version_ = std::stoi(val);
			else if(key == "infile") infile_ = val;
			else if(key == "detfile") detfile_ = val;
			else if(key == "treefile") treefile_ = val;
			else if(key == "histofile") histofile_ = val;
		}
		else if(section == "TargetLosses"){
			if(key.rfind("targetLoss_par",0) == 0){
				int idx = key.back() - '0';
				if(idx >= 1 && idx <= 4) targetLoss_par_[idx-1] = val;
			}
		}
		else if(section == "DeadLayerLosses"){
			if(key.rfind("deadLayerLoss_par",0) == 0){
				int idx = key.back() - '0';
				if(idx >= 1 && idx <= 4) deadLayerLoss_par_[idx-1] = val;
			}
		}
		else if(section == "Beamspot"){
			if(key == "profile") beam_profile_ = val;
			else if(key == "parX") beam_parX_ = std::stod(val);
			else if(key == "parY") beam_parY_ = std::stod(val);
			else if(key == "beam_offsetX") beam_offsetX_ = std::stod(val);
			else if(key == "beam_offsetY") beam_offsetY_ = std::stod(val);
		}
		else if(section == "MetaData"){
			if(key == "reaction") reaction_ = val;
			else if(key == "beam_energy") beam_energy_ = std::stod(val);
			else if(key == "recoil_excitation_energy") recoil_excitation_energy_ = std::stod(val);
		}
	}


	infile.close();
	return true;

}

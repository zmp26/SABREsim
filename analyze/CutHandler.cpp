#include "CutHandler.h"
#include <limits>
#include <cmath>

CutHandler::~CutHandler(){
	for(auto const& [key, cut] : fCuts2D) if(cut) delete cut;

	fCuts2D.clear();
	fCuts1D.clear();
}

void CutHandler::AddCut1D(const TString& name, double min, double max, bool isThreshold, bool isUpperLimit){
	Cut1D cut;
	cut.isThreshold = isThreshold;
	cut.isUpperLimit = isUpperLimit;

	cut.minVal = -1.0e30;
	cut.maxVal = 1.0e30;

	if(isThreshold) cut.minVal = min;

	if(isUpperLimit) cut.maxVal = max;

	if(!isThreshold && !isUpperLimit){
		cut.minVal = min;
		cut.maxVal = max;
	}

	fCuts1D[name] = cut;
}

//LoadCut2D() loads into fCuts2D[analysisKey] a single TCutG object named cutName from file located at fileName
bool CutHandler::LoadCut2D(const TString& fileName, const TString& cutName, const TString& analysisKey){
	TFile* file = TFile::Open(fileName, "READ");
	if(!file || file->IsZombie()){
		std::cerr << "CutHandler::LoadCut2D Error: Could not open file " << fileName << std::endl;
		return false;
	}

	TCutG* filecut = (TCutG*)file->Get(cutName);
	if(!filecut){
		std::cerr << "CutHandler::LoadCut2D Error: Could not find TCutG " << cutName << " in " << fileName << std::endl;
		file->Close();
		return false;
	}

	TCutG* cut = (TCutG*)filecut->Clone();

	//cut->SetDirectory(0);
	//fCuts2D.push_back(cut);
	fCuts2D[analysisKey] = cut;

	file->Close();
	delete file;
	return true;
}

//LoadAllCuts2D() loads all cuts with names in cutNames that appear in the file at fileName
//keep this function commented out as I plan on saving different cuts to different files
//however, write another function that reads in all 1D and 2D cuts defined in a file? --> should be a list of relevant info or path to file with relevant info
// bool CutHandler::LoadAllCuts2D(const TString& fileName, const std::vector<TString>& cutNames){
// 	bool allLoaded = true;
// 	for(const auto& name : cutNames){
// 		if(!LoadCut2D(fileName, name)){
// 			allLoaded = false;
// 		}
// 	}
// 	return allLoaded;
// }

void CutHandler::AddDirectCut2D(TCutG* cut, const TString& analysisKey){
	if(cut){
		//cut->SetDirectory(0);
		fCuts2D[analysisKey] = (TCutG*)cut->Clone();
	}
}

bool CutHandler::CheckCut1D(const TString& name, double value) const {
	auto it = fCuts1D.find(name);
	if(it != fCuts1D.end()) return (value >= it->second.minVal && value <= it->second.maxVal);
	return true;
}

bool CutHandler::CheckCut2D(const TString& name, double x, double y) const {
	auto it = fCuts2D.find(name);
	if(it != fCuts2D.end() && it->second != nullptr){
		if(it->second->IsInside(x,y)){
			std::cout << "Point (" << x << ", " << y << ") is inside cut with name " << name << "\n";
			return true;
		} else {
			//std::cout << "Point (" << x << ", " << y << ") is not inside cut with name " << name << "\n";
			return false;
		}
	}
	return false;
}

bool CutHandler::PassAll1D(const std::map<TString, double>& values) const {
	for(const auto& [gateName, cut] : fCuts1D){
		auto it = values.find(gateName);
		if(it != values.end()){
			if(it->second < cut.minVal || it->second > cut.maxVal){
				return false;
			}
		}
	}
	return true;
}

bool CutHandler::PassAll2D(const std::map<TString, std::pair<double,double>>& eventPoints) const {
	for(const auto& [keyName, cut] : fCuts2D){
		if(!cut) continue;

		auto it = eventPoints.find(keyName);
		if(it != eventPoints.end()){
			if(!cut->IsInside(it->second.first, it->second.second)){
				//std::cout << "point (" << it->second.first << ", " << it->second.second << ") is not inside " << keyName << "!\n";
				return false;
			} else {
				//std::cout << "point (" << it->second.first << ", " << it->second.second << ") is inside " << keyName << "!\n";
			}
		}
	}
	//std::cout << "returning true...\n";
	return true;
}

bool CutHandler::LoadCutsFromConfig(const TString& configFilePath){
	std::ifstream infile(configFilePath.Data());
	if(!infile.is_open()){
		std::cerr << "CutHandler::LoadCutsFromConfig Error: Unable to open config file " << configFilePath << "\n";
		return false;
	}

	std::string line;
	bool allLoaded = true;
	while(std::getline(infile, line)){
		if(line.empty() || line[0] == '#') continue;

		std::stringstream ss(line);
		std::string mode;
		ss >> mode;

		if(mode == "1D"){
			std::string name;
			double min, max;
			bool isThresh, isUpper;
			if(ss >> name >> min >> max >> isThresh >> isUpper){
				AddCut1D(name.c_str(), min, max, isThresh, isUpper);
			} else {
				std::cerr << "CutHandler Parsing Error on 1D line: " << line << "\n";
				allLoaded = false;
			}
		} else if(mode == "2D"){
			std::string filePath, cutObjectName, analysisKeyName;
			if(ss >> filePath >> cutObjectName >> analysisKeyName){
				if(!LoadCut2D(filePath.c_str(), cutObjectName.c_str(), analysisKeyName.c_str())){
					allLoaded = false;
				}
			} else {
				std::cerr << "CutHandler Parsing Error on 2D line: " << line << "\n";
				allLoaded = false;
			}
		}
	}
	return allLoaded;
}

void CutHandler::PrintCuts() const {
	std::cout << "\n========================================\n";
	std::cout << "      CutHandler: Registered Gates      \n";
	std::cout << "========================================\n";
	std::cout << "---> 1D Gates:\n";

	for(const auto& [gateName, cut] : fCuts1D){
		if(cut.isThreshold && !cut.isUpperLimit){
			std::cout << " [Lower Threshold] " << gateName << " >= " << cut.minVal << "\n";
		} else if(!cut.isThreshold && cut.isUpperLimit){
			std::cout << " [Upper Limit]     " << gateName << " <= " << cut.maxVal << "\n";			
		} else {
			std::cout << " [Range]           " << gateName << " : [" << cut.minVal << ", " << cut.maxVal << "]\n";
		}
	}

	std::cout << "\n---> 2D Graphical Polygons:\n";
	for(const auto& [keyName, cut] : fCuts2D){
		if(cut){
			std::cout << " [TCutG]           " << keyName << " (" << cut->GetN() << " vertices)\n";
		}
	}
	std::cout << "========================================\n\n";
}

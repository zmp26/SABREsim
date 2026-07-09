#include "CutHandler.h"
#include <limits>
#include <cmath>

CutHandler::~CutHandler(){
	for(auto cut : fCuts2D){
		if(cut){
			delete cut;
		}
	}
	fCuts2D.clear();
	fCuts1D.clear();
}

void CutHandler::AddCut1D(const TString& name, double min, double max, bool isThreshold, bool isUpperLimit){
	Cut1D cut;
	cut.name = name;
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

	fCuts1D.push_back(cut);
}

bool CutHandler::LoadCut2D(const TString& fileName, const TString& cutName){
	TFile* file = TFile::Open(fileName, "READ");
	if(!file || file->IsZombie()){
		std::cerr << "CutHandler::LoadCut2D Error: Could not open file " << fileName << std::endl;
		return false;
	}

	TCutG* cut = (TCutG*)file->Get(cutName);
	if(!cut){
		std::cerr << "CutHandler::LoadCut2D Error: Could not find TCutG " << cutName << " in " << fileName << std::endl;
		file->Close();
		return false;
	}

	//cut->SetDirectory(0);
	fCuts2D.push_back(cut);

	file->Close();
	return true;
}

bool CutHandler::LoadAllCuts2D(const TString& fileName, const std::vector<TString>& cutNames){
	bool allLoaded = true;
	for(const auto& name : cutNames){
		if(!LoadCut2D(fileName, name)){
			allLoaded = false;
		}
	}
	return allLoaded;
}

void CutHandler::AddDirectCut2D(TCutG* cut){
	if(cut){
		//cut->SetDirectory(0);
		fCuts2D.push_back(cut);
	}
}

bool CutHandler::CheckCut1D(const TString& name, double value) const {
	for(const auto& cut : fCuts1D){
		if(cut.name == name){
			//boundaries pre-padded to +/- 1e30 if a flag is not set
			return (value >= cut.minVal && value <= cut.maxVal);
		}
	}
	return true;
}

bool CutHandler::CheckCut2D(const TString& name, double x, double y) const {
	for(const auto& cut : fCuts2D){
		if(cut && cut->GetName() == name){
			return cut->IsInside(x,y);
		}
	}
	return false;
}

bool CutHandler::PassAll1D(const std::map<TString, double>& values) const {
	for(const auto& cut : fCuts1D){
		auto it = values.find(cut.name);
		if(it != values.end()){
			if(it->second < cut.minVal || it->second > cut.maxVal){
				return false;
			}
		}
	}
	return true;
}

bool CutHandler::PassAll2D(const std::map<TString, std::pair<double,double>>& eventPoints) const {
	for(auto cut : fCuts2D){
		if(!cut) continue;

		auto it = eventPoints.find(cut->GetName());
		if(it != eventPoints.end()){
			double x = it->second.first;
			double y = it->second.second;

			if(!cut->IsInside(x,y)){
				return false;
			}
		}
	}
	return true;
}

void CutHandler::PrintCuts() const {
	std::cout << "\n========================================\n";
	std::cout << "      CutHandler: Registered Gates      \n";
	std::cout << "========================================\n";
	std::cout << "---> 1D Gates:\n";

	for(const auto& cut : fCuts1D){
		if(cut.isThreshold && !cut.isUpperLimit){
			std::cout << " [Lower Threshold] " << cut.name << " >= " << cut.minVal << "\n";
		} else if(!cut.isThreshold && cut.isUpperLimit){
			std::cout << " [Upper Limit]     " << cut.name << " <= " << cut.maxVal << "\n";			
		} else {
			std::cout << " [Range]           " << cut.name << " : [" << cut.minVal << ", " << cut.maxVal << "]\n";
		}
	}

	std::cout << "\n---> 2D Graphical Polygons:\n";
	for(const auto& cut : fCuts2D){
		if(cut){
			std::cout << " [TCutG]           " << cut->GetName() << " (" << cut->GetN() << " vertices)\n";
		}
	}
	std::cout << "========================================\n\n";
}

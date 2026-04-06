#include "SABRE_AngleMap.h"

#include <fstream>
#include <sstream>
#include <iostream>

SABRE_AngleMap::SABRE_AngleMap(const std::string& filename){
	LoadMapFromFile(filename);
}

SABRE_AngleMap::SABRE_AngleMap(const std::vector<std::string>& filenames){
	for(const auto &fn : filenames){
		LoadMapFromFile(fn);
		//std::cout << "anglemap size = " << fAngleMap.size() << std::endl;
	};
}

bool SABRE_AngleMap::LoadMapFromFile(const std::string& filename){
	std::ifstream infile(filename);
	if(!infile.is_open()){
		std::cerr << "SABRE_AngleMap: Failed to open file: " << filename << std::endl;
		return false;
	}

	std::string line;
	int ring, wedge;
	double theta, phi;

	while(std::getline(infile,line)){

		if(line.empty() || line[0] == '#' || line[0] == 'r'){
			continue;
		}

		if(ParseLine(line,ring,wedge,theta,phi)){
			ChannelKey key = std::make_pair(ring, wedge);
			AngleValue val = std::make_pair(theta, phi);
			//std::cout << "Loaded (ring, wedge) = (" << ring << ", " << wedge << ") with (theta, phi) = (" << theta << ", " << phi << ")\n";

			fAngleMap[key] = val;
		} else {
			std::cerr << "SABRE_AngleMap: Failed to parse line: " << line << std::endl;
		}

	}

	return true;
}

std::optional<std::pair<double,double>> SABRE_AngleMap::GetDetectorThetaPhi(int ringchan, int wedgechan) const {
	ChannelKey key = std::make_pair(ringchan,wedgechan);

	auto it = fAngleMap.find(key);
	if(it != fAngleMap.end()){
		return it->second;
	}

	return std::nullopt;
}

bool SABRE_AngleMap::HasEntry(int ringchan, int wedgechan) const {
	ChannelKey key = std::make_pair(ringchan, wedgechan);
	return fAngleMap.find(key) != fAngleMap.end();
}

void SABRE_AngleMap::Clear(){
	fAngleMap.clear();
}

bool SABRE_AngleMap::ParseLine(const std::string& line, int& ring, int& wedge, double& theta, double& phi) const {
	std::istringstream iss(line);

	if(!(iss >> ring >> wedge >> theta >> phi)){
		return false;
	}

	return true;
}
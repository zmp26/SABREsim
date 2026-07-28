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
	double theta, sigmatheta, phi, sigmaphi;

	while(std::getline(infile,line)){

		if(line.empty() || line[0] == '#' || line[0] == 'r'){
			continue;
		}

		if(ParseLine(line,ring,wedge,theta,sigmatheta,phi,sigmaphi)){
			ChannelKey key = std::make_pair(ring, wedge);

			std::pair<double,double> thetapair = std::make_pair(theta, sigmatheta);
			std::pair<double,double> phipair = std::make_pair(phi, sigmaphi);
			AngleValue val = std::make_pair(thetapair, phipair);

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

	auto it = fAngleMap.find(key); //it is of same type as entry in fAngleMap which is std::pair<std::pair<int,int>, std::pair<std::pair<double,double>,std::pair<double,double>>>
									//this seems complicated, but looks like this in practice: {{ring, wedge}, {{theta, sigmatheta}, {phi, sigmaphi}}};
									//so, you find the entry based on the unique (ring, wedge) pair and get the theta, sigmatheta, phi, and sigmaphi out
									//TODO: rewrite this to take a struct or something that makes it much more readable than a bunch of 'first's and 'second's
	if(it != fAngleMap.end()){
		std::pair<double,double> retpair = std::make_pair(it->second.first.first, it->second.second.first);
		return retpair;
	}

	return std::nullopt;
}

std::optional<std::pair<double,double>> SABRE_AngleMap::GetDetectorThetaPhiSigmas(int ringchan, int wedgechan) const {
	ChannelKey key = std::make_pair(ringchan, wedgechan);

	auto it = fAngleMap.find(key);//it is of same type as entry in fAngleMap which is std::pair<std::pair<int,int>, std::pair<std::pair<double,double>,std::pair<double,double>>>
									//this seems complicated, but looks like this in practice: {{ring, wedge}, {{theta, sigmatheta}, {phi, sigmaphi}}};
									//so, you find the entry based on the unique (ring, wedge) pair and get the theta, sigmatheta, phi, and sigmaphi out
									//TODO: rewrite this to take a struct or something that makes it much more readable than a bunch of 'first's and 'second's

	if(it != fAngleMap.end()){
		std::pair<double,double> retpair = std::make_pair(it->second.first.second, it->second.second.second);
		return retpair;
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

bool SABRE_AngleMap::ParseLine(const std::string& line, int& ring, int& wedge, double& theta, double& sigmatheta, double& phi, double& sigmaphi) const {
	std::istringstream iss(line);

	if(!(iss >> ring >> wedge >> theta >> sigmatheta >> phi >> sigmaphi)){
		return false;
	}

	return true;
}
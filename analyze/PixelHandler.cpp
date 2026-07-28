#include "PixelHandler.h"

#include <fstream>
#include <sstream>
#include <iostream>

PixelHandler::PixelHandler(const std::string& filename){
	LoadFromFile(filename);
}

PixelHandler::PixelHandler(const std::vector<std::string>& filenames){
	for(const auto &fn : filenames){
		LoadFromFile(fn);
	}
}

bool PixelHandler::LoadFromFile(const std::string& filename){
	std::ifstream infile(filename);
	if(!infile.is_open()){
		std::cerr << "PixelHandler: Failed to open file: " << filename << std::endl;
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

			Pixel p;
			p.ring = ring;
			p.wedge = wedge;
			p.theta = theta;
			p.sigmatheta = sigmatheta;
			p.phi = phi;
			p.sigmaphi = sigmaphi;
			p.sigmaE = 0.05;//default 50 keV resolution for all pixels until determined later

			fAngleMap[key] = p;
		} else {
			std::cerr << "PixelHandler: Failed to parse line: " << line << std::endl;
		}

	}

	return true;
}

std::optional<std::pair<double,double>> PixelHandler::GetDetectorThetaPhi(int ringchan, int wedgechan) const {
	ChannelKey key = std::make_pair(ringchan,wedgechan);

	auto it = fAngleMap.find(key);
	if(it != fAngleMap.end()){
		std::pair<double,double> retpair = std::make_pair(it->second.theta, it->second.phi);
		return retpair;
	}

	return std::nullopt;
}

std::optional<std::pair<double,double>> PixelHandler::GetDetectorThetaPhiSigmas(int ringchan, int wedgechan) const {
	ChannelKey key = std::make_pair(ringchan, wedgechan);

	auto it = fAngleMap.find(key);
	if(it != fAngleMap.end()){
		std::pair<double,double> retpair = std::make_pair(it->second.sigmatheta, it->second.sigmaphi);
		return retpair;
	}

	return std::nullopt;
}

std::optional<double> PixelHandler::GetDetectorSigmaE(int ringchan, int wedgechan) const {
	ChannelKey key = std::make_pair(ringchan, wedgechan);

	auto it = fAngleMap.find(key);
	if(it != fAngleMap.end()){
		return it->second.sigmaE;
	}

	return std::nullopt;
}

std::optional<Pixel*> PixelHandler::GetDetectorPixel(int ringchan, int wedgechan) {
	ChannelKey key = std::make_pair(ringchan, wedgechan);

	auto it = fAngleMap.find(key);
	if(it != fAngleMap.end()){
		return &(it->second);
	}

	return std::nullopt;
}

bool PixelHandler::HasEntry(int ringchan, int wedgechan) const {
	ChannelKey key = std::make_pair(ringchan, wedgechan);
	return fAngleMap.find(key) != fAngleMap.end();
}

void PixelHandler::Clear(){
	fAngleMap.clear();
}

bool PixelHandler::ParseLine(const std::string& line, int& ring, int& wedge, double& theta, double& sigmatheta, double& phi, double& sigmaphi) const {
	std::istringstream iss(line);

	if(!(iss >> ring >> wedge >> theta >> sigmatheta >> phi >> sigmaphi)){
		return false;
	}

	return true;
}
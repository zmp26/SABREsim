#include "SABRE_ChannelMap.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

SABRE_ChannelMap::SABRE_ChannelMap(const std::string& filename){
	LoadMap(filename);
}

void SABRE_ChannelMap::LoadMap(const std::string& filename){
	//reads in from config/ChannelMap_Feb2021_SABRE.txt
	std::ifstream file(filename);
	if(!file.is_open()){
		throw std::runtime_error("Could not open channel map file: " + filename);
	}

	std::string line;
	int maxChannel = -1;

	std::string discard;
	//skip first two lines
	std::getline(file, discard);
	std::getline(file, discard);

	//first pass to extract maxChannel for memory purposes
	while(std::getline(file,line)){
		if(line.empty()) continue;
		std::istringstream ss(line);
		int globalChannel;
		ss >> globalChannel;
		if(globalChannel > maxChannel) maxChannel = globalChannel;
	}

	channelMap_.resize(maxChannel+1);

	file.clear();
	file.seekg(0,std::ios::beg);

	//skip first two lines
	std::getline(file, discard);
	std::getline(file, discard);

	//second pass to extract data
	while(std::getline(file,line)){
		if(line.empty()) continue;

		std::istringstream ss(line);
		int globalChannel;
		int detectorID;
		std::string detectorType;
		std::string detectorPart;

		ss >> globalChannel >> detectorID >> detectorType >> detectorPart;

		if(globalChannel >= 0 && globalChannel < static_cast<int>(channelMap_.size())){
			ChannelInfo info;
			info.detectorID = detectorID;
			info.detectorType = detectorType;
			info.detectorPart = detectorPart;
			info.valid = (detectorID != -1);

			channelMap_[globalChannel] = info;
		}
	}

	file.close();

}

const ChannelInfo& SABRE_ChannelMap::GetChannelInfo(int globalChannel) const {
	if(globalChannel < 0 || globalChannel >= static_cast<int>(channelMap_.size())){
		static ChannelInfo invalid = {-1, "INVALID", "INVALID", false};
		return invalid;
	}

	return channelMap_[globalChannel];
}
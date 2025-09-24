#include "SABRE_ChannelMap.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

SABRE_ChannelMap::SABRE_ChannelMap(const std::string& filename){
	std::ifstream file(filename);
	if(!file){
		throw std::runtime_error("Failed to open file: " + filename);
	}

	std::string line;
	int lineNumber = 0;

	while(lineNumber < 3 && std::getline(file, line)){
		lineNumber += 1;
	}

	while(std::getline(file, line)){
		if(line.empty()) continue;

		std::istringstream iss(line);
		int global_channel, detectorID_number;
		std::string detectorType_identifier, detectorPart_identifier;

		if(!(iss >> global_channel >> detectorID_number >> detectorType_identifier >> detectorPart_identifier)){
			throw std::runtime_error("Failed to parse line: " + line);
		}

		SABRE_ChannelMap::DetectorKey key = std::make_tuple(detectorID_number, detectorType_identifier, detectorPart_identifier);
		keyToGlobalChannel[key] = global_channel;
		globalChannelToKey[global_channel] = key;
	}
}

int SABRE_ChannelMap::getGlobalChannel(int detectorID, const std::string& type, const std::string& part) const {
	SABRE_ChannelMap::DetectorKey key = std::make_tuple(detectorID, type, part);
	auto it = keyToGlobalChannel.find(key);
	if(it == keyToGlobalChannel.end()){
		throw std::runtime_error("Detector key not found");
	}
	return it->second;
}

SABRE_ChannelMap::DetectorKey SABRE_ChannelMap::getDetectorKey(int global_channel) const {
	auto it = globalChannelToKey.find(global_channel);
	if(it == globalChannelToKey.end()){
		throw std::runtime_error("Global channel not found");
	}
	return it->second;
}
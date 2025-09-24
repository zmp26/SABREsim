#ifndef SABRE_CHANNEL_MAP_H
#define SABRE_CHANNEL_MAP_H

#include <string>
#include <unordered_map>
#include <tuple>

class SABRE_ChannelMap {
public:
	using DetectorKey = std::tuple<int, std::string, std::string>;

	SABRE_ChannelMap(const std::string& filename);

	int getGlobalChannel(int detectorID, const std::string& type, const std::string& part) const;
	DetectorKey getDetectorKey(int global_channel) const;

private:
	std::unordered_map<DetectorKey, int> keyToGlobalChannel;
	std::unordered_map<int, DetectorKey> globalChannelToKey;
};

namespace std {
	template<>
	struct hash<SABRE_ChannelMap::DetectorKey> {
		size_t operator()(const SABRE_ChannelMap::DetectorKey& key) const {
			auto h1 = std::hash<int>{}(std::get<0>(key));
			auto h2 = std::hash<std::string>{}(std::get<1>(key));
			auto h3 = std::hash<std::string>{}(std::get<2>(key));
			return h1 ^ (h2 << 1) ^ (h3 << 2); //bit shifting to help spread out hash values -> reduce risk of collisions when similar values get used!
		}
	};
}


#endif
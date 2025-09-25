#ifndef SABRE_CHANNEL_MAP_H
#define SABRE_CHANNEL_MAP_H

#include <string>
#include <vector>
#include <unordered_map>

struct ChannelInfo {
	int detectorID;
	std::string detectorType;
	std::string detectorPart;
	bool valid;
};

class SABRE_ChannelMap{
public:
	explicit SABRE_ChannelMap(const std::string& filename);

	const ChannelInfo& GetChannelInfo(int globalChannel) const;

private:
	std::vector<ChannelInfo> channelMap_;
	void LoadMap(const std::string& filename);
};

#endif
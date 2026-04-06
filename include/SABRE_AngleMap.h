#ifndef SABRE_ANGLEMAP_H
#define SABRE_ANGLEMAP_H

#include <map>
#include <utility>
#include <string>
#include <optional>
#include <vector>

class SABRE_AngleMap{
public:
	SABRE_AngleMap() = default;
	explicit SABRE_AngleMap(const std::string& filename);
	explicit SABRE_AngleMap(const std::vector<std::string>& filenames);

	~SABRE_AngleMap() = default;

	bool LoadMapFromFile(const std::string& filename);

	std::optional<std::pair<double,double>> GetDetectorThetaPhi(int ringchan, int wedgechan) const;

	bool HasEntry(int ringchan, int wedgechan) const;
	void Clear();


private:
	using ChannelKey = std::pair<int,int>;
	using AngleValue = std::pair<double,double>;

	std::map<ChannelKey, AngleValue> fAngleMap;

	bool ParseLine(const std::string& line, int& ring, int& wedge, double& theta, double& phi) const;
};



#endif//SABRE_ANGLEMAP_H
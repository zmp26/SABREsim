#ifndef PIXELHANDLER_H
#define PIXELHANDLER_H

#include <map>
#include <utility>
#include <string>
#include <optional>
#include <vector>

struct Pixel {
	int ring, wedge;
	double theta, sigmatheta; //theta is mean value and sigmatheta is std devition from anglemaps/sigmas/rays.cpp
	double phi, sigmaphi; //phi is mean value and sigmaphi is std deviation from anglemaps/sigmas/rays.cpp
	double sigmaE; //sigmaE is energy resolution of pixel determined from SABRE singles
};

class PixelHandler{
public:
	PixelHandler() = default;
	explicit PixelHandler(const std::string& filename);
	explicit PixelHandler(const std::vector<std::string>& filenames);

	~PixelHandler() = default;

	bool LoadFromFile(const std::string& filename);

	std::optional<std::pair<double,double>> GetDetectorThetaPhi(int ringchan, int wedgechan) const;
	std::optional<std::pair<double,double>> GetDetectorThetaPhiSigmas(int ringchan, int wedgechan) const;
	std::optional<double> GetDetectorSigmaE(int ringchan, int wedgechan) const;

	std::optional<Pixel*> GetDetectorPixel(int ringchan, int wedgechan);

	bool HasEntry(int ringchan, int wedgechan) const;
	void Clear();

private:
	using ChannelKey = std::pair<int,int>;
	std::map<ChannelKey, Pixel> fAngleMap;

	bool ParseLine(const std::string& line, int& ring, int& wedge, double& theta, double& sigmatheta, double& phi, double& sigmaphi) const;
};


#endif //PIXELHANDLER_H
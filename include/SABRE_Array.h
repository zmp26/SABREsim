#ifndef SABRE_ARRAY_H
#define SABRE_ARRAY_H

#include <vector>
#include "SABRE_Detector.h"
#include <string>
#include <optional>


class SABRE_Array{
public:
	SABRE_Array();
	~SABRE_Array();

	void Clear() { SABRE_Array_.clear(); };

	void push_back(SABRE_Detector* det) { SABRE_Array_.push_back(det); };

	size_t size() { return SABRE_Array_.size(); };

	SABRE_Detector* at(size_t i){ return SABRE_Array_.at(i); };

	std::optional<std::pair<double,double>> GetDetectorThetaPhi(int ringchan, int wedgechan);




private:

	std::vector<SABRE_Detector*> SABRE_Array_;

	const std::vector<std::string> filenames = {
		"anglemaps/SABRE0_phi306_anglemap.txt",
		"anglemaps/SABRE1_phi18_anglemap.txt",
		"anglemaps/SABRE2_phi234_anglemap.txt",
		"anglemaps/SABRE3_phi162_anglemap.txt",
		"anglemaps/SABRE4_phi90_anglemap.txt"
	};

	SABRE_AngleMap *anglemap;

};





#endif//SABRE_ARRAY_H
#include "SABRE_Array.h"
#include "SABRE_AngleMap.h"

SABRE_Array::SABRE_Array(){
	anglemap = new SABRE_AngleMap(filenames);
}

SABRE_Array::~SABRE_Array(){};

std::optional<std::pair<double,double>> SABRE_Array::GetDetectorThetaPhi(int ringchan, int wedgechan){
	return anglemap->GetDetectorThetaPhi(ringchan, wedgechan);
}
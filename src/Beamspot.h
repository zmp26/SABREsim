//code taken from "detector_stm_v2.h" provided to me by STM
//renaming class from TBeamspot to just Beamspot:

#ifndef BEAMSPOT_H
#define BEAMSPOT_H

using namespace std;

class Beamspot
{
private:
	float fXMax;
	float fYMax;
	float fDx;
	float fDy;
	float fXOffset;
	float fYOffset;
	float fTheta;
	float fPhi;
	const double DEGRAD = 0.017453293;
public:
	void SetXMax(float xmax) { fXMax = xmax; }
	void SetYMax(float ymax) { fYMax = ymax; }
	void SetXOffset(float xoffset) { fXOffset = xoffset; }
	void SetYOffset(float yoffset) { fYOffset = yoffset; }
	void Spread(void);
	void Set(float z, float theta, float phi);
	float GetTheta(void) { return fTheta; }
	float GetPhi(void) { return fPhi; }
	float GetDx(void) { return fDx; }
	float GetDy(void) { return fDy; }
	void Print(void);
};

#endif
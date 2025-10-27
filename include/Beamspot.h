//code taken from "detector_stm_v2.h" provided to me by STM
//renaming class from TBeamspot to just Beamspot:

#ifndef BEAMSPOT_H
#define BEAMSPOT_H

// using namespace std;

// #include "Vec3.h"

// class Beamspot
// {
// private:
// 	double fXMax;
// 	double fYMax;
// 	double fDx;
// 	double fDy;
// 	double fXOffset;
// 	double fYOffset;
// 	double fTheta;
// 	double fPhi;
// 	const double DEGRAD = 0.017453293;
// public:
// 	void SetXMax(double xmax) { fXMax = xmax; }
// 	void SetYMax(double ymax) { fYMax = ymax; }
// 	void SetXOffset(double xoffset) { fXOffset = xoffset; }
// 	void SetYOffset(double yoffset) { fYOffset = yoffset; }
// 	void SetTheta(double theta) { fTheta = theta; }
// 	void SetPhi(double phi) { fPhi = phi; }
// 	void Spread(void);
// 	void Set(double z, double theta, double phi);
// 	bool Set(const Vec3& detectorPoint, const Vec3& detectorNormal, double theta, double phi);
// 	double GetTheta(void) { return fTheta; }
// 	double GetPhi(void) { return fPhi; }
// 	double GetDx(void) { return fDx; }
// 	double GetDy(void) { return fDy; }
// 	void Print(void);
// };

#include "BeamProfile.h"
#include "Vec3.h"
#include "TString.h"

class Beamspot {
public:
	Beamspot();
	~Beamspot();

	//Disable copy via constructor explicitly deleting
	Beamspot(const Beamspot&) = delete;
	Beamspot& operator=(const Beamspot&) = delete;

	Beamspot(Beamspot&& other) noexcept;
	Beamspot& operator=(Beamspot&& other) noexcept;

	void SetProfile(BeamProfile* p);
	void SetBeamAxisOffset(double xoff, double yoff);

	Vec3 GeneratePoint(double z = 0.0) const;

	double GetXOffset(){ return xOffset; }
	double GetYOffset(){ return yOffset; }

private:
	BeamProfile* profile = nullptr;
	double xOffset = 0.;
	double yOffset = 0.;

};

#endif
//code taken from "detector_stm_v2.cpp" provided to me by STM
//renaming class from TBeamspot to just Beamspot:
//using namespace std;

// #include <iostream>
// #include <cmath>

// #include "Beamspot.h"
// #include "Vec3.h"

// void Beamspot::Spread(void){
// 	fDx = 2*fXMax*((rand() % 10000)/10000.) - fXMax + fXOffset;
// 	fDy = 2*fYMax*((rand() % 10000)/10000.) - fYMax + fYOffset;
// }

// void Beamspot::Set(double z, double theta, double phi){
// 	double r, r2, theta2, phi2, x, y, x2, y2;
// 	if((fXMax==0 || fYMax==0)&&(fXOffset==0 && fYOffset==0)){
// 		fTheta = theta;
// 		fPhi = phi;
// 		return;
// 	}
// 	r=z/std::cos(theta*DEGRAD);
// 	x=r*std::sin(theta*DEGRAD)*std::cos(phi*DEGRAD);
// 	y=r*std::sin(theta*DEGRAD)*std::sin(phi*DEGRAD);
// 	x2=x+fDx;
// 	y2=y+fDy;
// 	if(r>0){
// 		r2 = std::sqrt(x2*x2 + y2*y2 + z*z);
// 	} else {
// 		r2=-1*std::sqrt(x2*x2 + y2*y2 + z*z);
// 	}
// 	theta2=std::acos(z/r2)/DEGRAD;
// 	phi2=std::acos(y2,x2)/DEGRAD;
// 	fTheta=theta2;
// 	if(phi2>0){
// 		fPhi=phi2;
// 	} else {
// 		fPhi=phi2+360.;
// 	}
// }

// void Beamspot::Print(void){
// 	std::cout << "Interaction point is: Dx = " << fDx << " Dy = " << fDy << std::endl;
// }

// bool Beamspot::Set(const Vec3& detectorPoint, const Vec3& detectorNormal, double theta, double phi){

// 	double thetarad = theta*DEGRAD;
// 	double phirad = phi*DEGRAD;

// 	Vec3 trajectory(std::sin(thetarad)*std::cos(phirad), std::sin(thetarad)*std::sin(phirad), std::cos(thetarad));
// 	Vec3 origin(GetDx(), GetDy(), 0);

// 	Vec3 p0 = detectorPoint;

// 	double denom = detectorNormal.Dot(trajectory);
// 	if(fabs(denom) < 1e-6){
// 		//cerr << "Warning: trajectory direction is parallel to detector plane\n";
// 		fTheta = theta;
// 		fPhi = phi;
// 		return false;
// 	}

// 	double t = detectorNormal.Dot(detectorPoint - origin) / denom;//parameter t for when origin + trajectory*t intersects the plane!
// 	Vec3 intersection = origin + trajectory*t;//set the intersection point

// 	Vec3 newdir = intersection;
// 	double r = newdir.Mag();
// 	double newtheta = std::acos(newdir.GetZ()/r)/DEGRAD;
// 	double newphi = std::acos(newdir.GetY(), newdir.GetX())/DEGRAD;
// 	if(newphi < 0) newphi += 360.;

// 	//use original theta, phi if no offset:
// 	if((fXMax==0 || fYMax==0) && (fXOffset==0 && fYOffset==0)){
// 		fTheta = theta;
// 		fPhi = phi;
// 		//std::cout << "no offset found, using original theta and phi" << std::endl;
// 	} else {
// 		fTheta = newtheta;
// 		fPhi = newphi;
// 		//std::cout << "updating beamspot with new theta and phi" << std::endl;
// 	}

// 	return true;
// }

#include "Beamspot.h"

Beamspot::Beamspot() = default;

Beamspot::~Beamspot() {
	delete profile;
}

Beamspot::Beamspot(Beamspot&& other) noexcept {
	profile = other.profile;
	xOffset = other.xOffset;
	yOffset = other.yOffset;
	other.profile = nullptr;
}

Beamspot& Beamspot::operator=(Beamspot&& other) noexcept{
	if(this != &other){
		delete profile;
		profile = other.profile;
		xOffset = other.xOffset;
		yOffset = other.yOffset;
		other.profile = nullptr;
	}

	return *this;
}

void Beamspot::SetProfile(BeamProfile* p){
	delete profile;
	profile = p;
}

void Beamspot::SetBeamAxisOffset(double xoff, double yoff){
	xOffset = xoff;
	yOffset = yoff;
}

Vec3 Beamspot::GeneratePoint(double z) const {
	if(!profile){
		return Vec3(xOffset, yOffset, z);//fallback of beam axis if no profile
	}

	auto[dx,dy] = profile->Sample();
	return Vec3(dx+xOffset, dy+yOffset, z);
}
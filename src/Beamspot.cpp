//code taken from "detector_stm_v2.cpp" provided to me by STM
//renaming class from TBeamspot to just Beamspot:
using namespace std;

#include <iostream>
#include <cmath>

#include "Beamspot.h"
#include "Vec3.h"

void Beamspot::Spread(void){
	fDx = 2*fXMax*((rand() % 10000)/10000.) - fXMax + fXOffset;
	fDy = 2*fYMax*((rand() % 10000)/10000.) - fYMax + fYOffset;
}

void Beamspot::Set(double z, double theta, double phi){
	double r, r2, theta2, phi2, x, y, x2, y2;
	if((fXMax==0 || fYMax==0)&&(fXOffset==0 && fYOffset==0)){
		fTheta = theta;
		fPhi = phi;
		return;
	}
	r=z/cos(theta*DEGRAD);
	x=r*sin(theta*DEGRAD)*cos(phi*DEGRAD);
	y=r*sin(theta*DEGRAD)*sin(phi*DEGRAD);
	x2=x+fDx;
	y2=y+fDy;
	if(r>0){
		r2 = sqrt(x2*x2 + y2*y2 + z*z);
	} else {
		r2=-1*sqrt(x2*x2 + y2*y2 + z*z);
	}
	theta2=acos(z/r2)/DEGRAD;
	phi2=atan2(y2,x2)/DEGRAD;
	fTheta=theta2;
	if(phi2>0){
		fPhi=phi2;
	} else {
		fPhi=phi2+360.;
	}
}

void Beamspot::Print(void){
	cout << "Interaction point is: Dx = " << fDx << " Dy = " << fDy << endl;
}

bool Beamspot::Set(const Vec3& detectorPoint, const Vec3& detectorNormal, double theta, double phi){

	double thetarad = theta*DEGRAD;
	double phirad = phi*DEGRAD;

	Vec3 trajectory(sin(thetarad)*cos(phirad), sin(thetarad)*sin(phirad), cos(thetarad));
	Vec3 origin(GetDx(), GetDy(), 0);

	Vec3 p0 = detectorPoint;

	double denom = detectorNormal.Dot(trajectory);
	if(fabs(denom) < 1e-6){
		//cerr << "Warning: trajectory direction is parallel to detector plane\n";
		fTheta = theta;
		fPhi = phi;
		return false;
	}

	double t = detectorNormal.Dot(detectorPoint - origin) / denom;//parameter t for when origin + trajectory*t intersects the plane!
	Vec3 intersection = origin + trajectory*t;//set the intersection point

	Vec3 newdir = intersection;
	double r = newdir.Mag();
	double newtheta = acos(newdir.GetZ()/r)/DEGRAD;
	double newphi = atan2(newdir.GetY(), newdir.GetX())/DEGRAD;
	if(newphi < 0) newphi += 360.;

	//use original theta, phi if no offset:
	if((fXMax==0 || fYMax==0) && (fXOffset==0 && fYOffset==0)){
		fTheta = theta;
		fPhi = phi;
		//cout << "no offset found, using original theta and phi" << endl;
	} else {
		fTheta = newtheta;
		fPhi = newphi;
		//cout << "updating beamspot with new theta and phi" << endl;
	}

	return true;
}
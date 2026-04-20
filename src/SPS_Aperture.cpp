#include "SPS_Aperture.h"

SPS_Aperture::SPS_Aperture(double r,
						   double thetaMin, double thetaMax,
						   double phiMin, double phiMax,
						   double sigmaE, double sigmaTheta, double sigmaPhi)
	: radius(r), thetaMin(thetaMin), thetaMax(thetaMax), phiMin(phiMin), phiMax(phiMax),
	  sigmaE(sigmaE), sigmaTheta(sigmaTheta), sigmaPhi(sigmaPhi)
{
	std::random_device rd;
	gen.seed(rd());

	std::cout << "Theta range: [" << thetaMin << ", " << thetaMax << "]\n";
	std::cout << "Phi range: [" << phiMin << ", " << phiMax << "]\n";
}

SPS_Aperture::~SPS_Aperture() {}

/*
 *	IsDetected() solves for intersection of a ray with trajectory traj from some
 *  point origin with a sphere within some bounds [thetaMin, thetaMax] and 
 *  [phiMin, phiMax] (these are set in the class)
 *
 *	Returns true if intersection occurs, false otherwise
 *	
 *  Parameterizes as follows:
 *		Ray: 		P(t) = origin + t*unit_traj
 *		Sphere: 	|P(t)|^2 = radius^2
 *
 *	Solves for the t(s) in the ray that intersects a sphere centered at origin
 *  and determins if within bounds
 *
 *  Quadratic: 
 *		(O + td) * (O + td) = R^2 where O is the reaction origin and d the trajectory,
 *		and R the distance from target plane to aperture (* is dot product here)
 *		==>		(d^2)t^2 + (2(O*d))t + (O^2 - R^2) = 0
 *		==>		A = d^2 	B = 2(O*d)		C = (O^2 - R^2)
 *		
 *	Solved by:
 *		t = (-2(O*d) +/- SQRT[ 4(O*d) - 4d^2(O^2 - R^2) ]) / (2d^2)
 *		NOTE: Code chooses positive root
 */
bool SPS_Aperture::IsDetected(const  Vec3& traj, const Vec3& origin) const {

	Vec3 d = traj.Unit();

	double A = d.Mag2(); // = 1 because d is unit vector
	double B = 2. * origin.Dot(d);
	double C = origin.Mag2() - radius*radius;

	double disc = (B*B) - 4.*A*C;

	//if disc < 0, trajectory misses sphere entirely... this is unlikely and should never happen, but handle anyway
	if(disc < 0) return false;

	double t = (-B + std::sqrt(disc)) / (2.*A);
	if(t<0) return false;

	Vec3 intersectPoint = origin + (d*t);

	double pointTheta = intersectPoint.GetTheta()*RADDEG;
	double pointPhi = intersectPoint.GetPhi()*RADDEG;
	//if(pointPhi < 0.) pointPhi += 360.;

	// std::cout << "intersectPoint: X = " << intersectPoint.GetX() << ", Y = " << intersectPoint.GetY() << ", Z = " << intersectPoint.GetZ() << "\n" 
	// 		  << "			  Theta = " << pointTheta << ", Phi = " << pointPhi << "\n\n";

	if( pointTheta >= thetaMin && pointTheta <= thetaMax && ((pointPhi >= 0. && pointPhi < phiMax) || (pointPhi <= 0. && pointPhi > phiMin)) ){
		return true;
	}

	return false;
}

double SPS_Aperture::GetSmearedEnergy(double energy) const {
	return ApplyGaussian(energy, sigmaE);
}

double SPS_Aperture::GetSmearedTheta(double theta) const {
	return ApplyGaussian(theta, sigmaTheta);
}

double SPS_Aperture::GetSmearedPhi(double phi) const {
	return ApplyGaussian(phi, sigmaPhi);
}

double SPS_Aperture::ApplyGaussian(double value, double sigma) const {
	if(sigma <= 0.) return value;

	std::normal_distribution<double> dist(value,sigma);

	return dist(gen);
}
/*
	This class takes heavy inspiration from "SABRE_Detector" written originally by KGH at FSU and then rewritten by GWM at FSU.

	This class represents a single Micron MMM detector in the SABRE Array at the FSU accelerator facility.
*/
/*

  Class which represents a single MMM detector in the SABRE array at FSU. Origial code by KGH, re-written by
  GWM.

  Distances in meters, angles in radians.

  The channel arrays have four points, one for each corner. The corners are
  as follows, as if looking BACK along beam (i.e. from the target's pov):

  0---------------------1
  |                     |
  |                     |      x
  |                     |      <-----
  |                     |      		|
  |                     |      		|
  3---------------------2      		y
                               (z is hence positive along beam direction) 

  The channel numbers, also as looking back from target pov, are:

  >> rings are 0 -- 15 from inner to outer:

    15 -------------------
    14 -------------------
    13 -------------------
       .
       .
       .
     2 -------------------
     1 -------------------
     0 -------------------

  >> wedges are 0 -- 7 moving counterclockwise:

      7 6 ... 1 0
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |


  >> Note that the detector starts centered on the x-axis (central phi = 0) untilted, and then is rotated to wherever the frick
  	 it is supposed to go; phi = 90 is centered on y axis, pointing down towards the bottom of the scattering chamber

  -- gwm, Dec 2020; based on the og code from kgh

*/
#include <random>
#include <cmath>
#include <fstream>
#include <iostream>
#include "SABRE_Detector.h"
#include "Vec3.h"

SABRE_Detector::SABRE_Detector() :
m_Router(0.1351), m_Rinner(0.0326), m_deltaPhi_flat(54.4*deg2rad), m_phiCentral(0.0), m_tilt(0.0), m_translation(0.,0.,0.), m_norm_flat(0,0,1.0)
{
	m_YRot.SetAngle(m_tilt);
	m_ZRot.SetAngle(m_phiCentral);

	//Initialize the coordinate arrays
	m_ringCoords_flat.resize(m_nRings);
	m_ringCoords_tilt.resize(m_nRings);
	m_wedgeCoords_flat.resize(m_nWedges);
	m_wedgeCoords_tilt.resize(m_nWedges);
	for(int i=0; i<m_nRings; i++) {
		m_ringCoords_flat[i].resize(4);
		m_ringCoords_tilt[i].resize(4);
	}
	for(int i=0; i<m_nWedges; i++) {
		m_wedgeCoords_flat[i].resize(4);
		m_wedgeCoords_tilt[i].resize(4);
	}

	m_deltaR_flat = m_Router - m_Rinner;
	m_deltaR_flat_ring = m_deltaR_flat/m_nRings;

	CalculateCorners();
}

SABRE_Detector::SABRE_Detector(double Rin, double Rout, double deltaPhi_flat, double phiCentral, double tiltFromVert, double zdist, double xdist, double ydist) :
m_Router(Rout), m_Rinner(Rin), m_deltaPhi_flat(deltaPhi_flat), m_phiCentral(phiCentral), m_tilt(tiltFromVert), m_translation(xdist, ydist, zdist), m_norm_flat(0,0,1.0)
{
	m_YRot.SetAngle(m_tilt);
	m_ZRot.SetAngle(m_phiCentral);

	//Initialize coordinate arrays
	m_ringCoords_flat.resize(m_nRings);
	m_ringCoords_tilt.resize(m_nRings);
	m_wedgeCoords_flat.resize(m_nWedges);
	m_wedgeCoords_tilt.resize(m_nWedges);
	for(int i=0; i<m_nRings; i++) {
		m_ringCoords_flat[i].resize(4);
		m_ringCoords_tilt[i].resize(4);
	}
	for(int i=0; i<m_nWedges; i++) {
		m_wedgeCoords_flat[i].resize(4);
		m_wedgeCoords_tilt[i].resize(4);
	}

	m_deltaR_flat = m_Router - m_Rinner;
	m_deltaR_flat_ring = m_deltaR_flat/m_nRings;
	m_deltaPhi_flat_wedge = m_deltaPhi_flat/m_nWedges;

	CalculateCorners();
}

SABRE_Detector::~SABRE_Detector() {}

void SABRE_Detector::CalculateCorners() {

	double x0, x1, x2, x3;
	double y0, y1, y2, y3;
	double z0, z1, z2, z3;

	//Generate flat ring corner coordinates
	for(int i=0; i<m_nRings; i++) {
		x0 = (m_Rinner + m_deltaR_flat_ring*(i+1))*std::cos(-m_deltaPhi_flat/2.0);
		y0 = (m_Rinner + m_deltaR_flat_ring*(i+1))*std::sin(-m_deltaPhi_flat/2.0);
		z0 = 0.0;
		m_ringCoords_flat[i][0].SetVectorCartesian(x0, y0, z0);

		x1 = (m_Rinner + m_deltaR_flat_ring*(i))*std::cos(-m_deltaPhi_flat/2.0);
		y1 = (m_Rinner + m_deltaR_flat_ring*(i))*std::sin(-m_deltaPhi_flat/2.0);
		z1 = 0.0;
		m_ringCoords_flat[i][1].SetVectorCartesian(x1, y1, z1);

		x2 = (m_Rinner + m_deltaR_flat_ring*(i))*std::cos(m_deltaPhi_flat/2.0);
		y2 = (m_Rinner + m_deltaR_flat_ring*(i))*std::sin(m_deltaPhi_flat/2.0);
		z2 = 0.0;
		m_ringCoords_flat[i][2].SetVectorCartesian(x2, y2, z2);

		x3 = (m_Rinner + m_deltaR_flat_ring*(i+1))*std::cos(m_deltaPhi_flat/2.0);
		y3 = (m_Rinner + m_deltaR_flat_ring*(i+1))*std::sin(m_deltaPhi_flat/2.0);
		z3 = 0.0;
		m_ringCoords_flat[i][3].SetVectorCartesian(x3, y3, z3);
	}

	//Generate flat wedge corner coordinates
	for(int i=0; i<m_nWedges; i++) {
		x0 = m_Router * std::cos(-m_deltaPhi_flat/2.0 + i*m_deltaPhi_flat_wedge);
		y0 = m_Router * std::sin(-m_deltaPhi_flat/2.0 + i*m_deltaPhi_flat_wedge);
		z0 = 0.0;
		m_wedgeCoords_flat[i][0].SetVectorCartesian(x0, y0, z0);

		x1 = m_Rinner * std::cos(-m_deltaPhi_flat/2.0 + i*m_deltaPhi_flat_wedge);
		y1 = m_Rinner * std::sin(-m_deltaPhi_flat/2.0 + i*m_deltaPhi_flat_wedge);
		z1 = 0.0;
		m_wedgeCoords_flat[i][1].SetVectorCartesian(x1, y1, z1);

		x2 = m_Rinner * std::cos(-m_deltaPhi_flat/2.0 + (i+1)*m_deltaPhi_flat_wedge);
		y2 = m_Rinner * std::sin(-m_deltaPhi_flat/2.0 + (i+1)*m_deltaPhi_flat_wedge);
		z2 = 0.0;
		m_wedgeCoords_flat[i][2].SetVectorCartesian(x2, y2, z2);

		x3 = m_Router * std::cos(-m_deltaPhi_flat/2.0 + (i+1)*m_deltaPhi_flat_wedge);
		y3 = m_Router * std::sin(-m_deltaPhi_flat/2.0 + (i+1)*m_deltaPhi_flat_wedge);
		z3 = 0.0;
		m_wedgeCoords_flat[i][3].SetVectorCartesian(x3, y3, z3);
	}

	//Generate tilted rings
	for(int i=0; i<m_nRings; i++) {
		for(int j=0; j<4; j++) {
			m_ringCoords_tilt[i][j] = TransformToTiltedFrame(m_ringCoords_flat[i][j]);
		}
	}

	//Generate tilted wedges
	for(int i=0; i<m_nWedges; i++) {
		for(int j=0; j<4; j++) {
			m_wedgeCoords_tilt[i][j] = TransformToTiltedFrame(m_wedgeCoords_flat[i][j]);
		}
	}
}

/*
	Given a unit vector (R=1, theta, phi) which corresponds to some particle's trajectory,
	determine whether that particle will intersect with this SABRE detector. If it does calculate
	the coordinates of the hit. The equation is as follows:
	
	Rz(eta)*Ry(psi)*Flat_vector(R', theta'=PI/2, phi') + translation = Tilted_vector(R, theta, phi)

	Where Rz is the ZRotation, Ry is the YRotation, F_vector is the vector of the coordinates in the flat detector frame,
	and Tilted_vector is the vector of the hit coordinates in the tilted frame. The theta and phi of the the Tilted_vector correspond
	to the input arguments of the function.

	It checks to deterime whether or not the particle hits within the borders (read: edges) of the SABRE detector, and does not account for
	the spacing between rings and wedges.

	!NOTE: This currently only applies to a configuration where there is no translation in x & y. The math becomes significantly messier in these cases.
	Also, don't use tan(). It's behavior near PI/2 makes it basically useless for these.
*/
Vec3 SABRE_Detector::GetTrajectoryCoordinates(double theta, double phi) {
	if(m_translation.GetX() != 0.0 || m_translation.GetY() != 0.0) {
		return Vec3();
	}

	//Calculate the *potential* phi in the flat detector
	double phi_numerator = std::cos(m_tilt)*(std::sin(phi)*std::cos(m_phiCentral) - std::sin(m_phiCentral)*std::cos(phi));
	double phi_denominator = std::cos(m_phiCentral)*std::cos(phi) + std::sin(m_phiCentral)*std::sin(phi);
	double phi_flat = std::atan2(phi_numerator, phi_denominator);
	if(phi_flat < 0) phi_flat += M_PI*2.0;

	//Calculate the *potential* R in the flat detector
	double r_numerator = m_translation.GetZ()*std::cos(phi)*std::sin(theta);
	double r_denominator = std::cos(phi_flat)*std::cos(m_phiCentral)*std::cos(m_tilt)*std::cos(theta) - std::sin(phi_flat)*std::sin(m_phiCentral)*std::cos(theta) - std::cos(phi_flat)*std::sin(m_tilt)*std::cos(phi)*std::sin(theta);
	double r_flat = r_numerator/r_denominator;

	//Calculate the distance from the origin to the hit on the detector
	double R_to_detector = (r_flat*std::cos(phi_flat)*std::sin(m_tilt) + m_translation.GetZ())/std::cos(theta);
	double xhit = R_to_detector*std::sin(theta)*std::cos(phi);
	double yhit = R_to_detector*std::sin(theta)*std::sin(phi);
	double zhit = R_to_detector*std::cos(theta);


	//Check to see if our flat coords fall inside the flat detector
	if(IsInside(r_flat, phi_flat)) {
		return Vec3(xhit, yhit, zhit);
	} else {
		return Vec3();
	}
}

/*
	Given a unit vector (R=1, theta, phi) which corresponds to some particle's trajectory,
	determine whether that particle will intersect with this SABRE detector. If it does determine 
	which ring and wedge the hit occurs in. The equation is as follows:
	
	Rz(eta)*Rx(psi)*Flat_vector(R', theta'=PI/2, phi') + translation = Tilted_vector(R, theta, phi)

	Where Rz is the ZRotation, Rx is the XRotation, F_vector is the vector of the coordinates in the flat detector frame,
	and Tilted_vector is the vector of the hit coordinates in the tilted frame. The theta and phi of the the Tilted_vector correspond
	to the input arguments of the function.

	Then using the flat coordinate R' and phi' determine which ring/wedge channels are hit. For precision purposes, the channel is not calculated, but
	rather found using comparisions. This method accounts for the spacing between rings and wedges.

	!NOTE: This currently only applies to a configuration where there is no translation in x & y. The math becomes significantly messier in these cases.
	Also, don't use tan(). It's behavior near PI/2 makes it basically useless for these.
*/
std::pair<int, int> SABRE_Detector::GetTrajectoryRingWedge(double theta, double phi) {
	if(m_translation.GetX() != 0.0 || m_translation.GetY() != 0.0) {
		return std::make_pair(-1, -1);
	}

	//Calculate the *potential* phi in the flat detector
	double phi_numerator = std::cos(m_tilt)*(std::sin(phi)*std::cos(m_phiCentral) - std::sin(m_phiCentral)*std::cos(phi));
	double phi_denominator = std::cos(m_phiCentral)*std::cos(phi) + std::sin(m_phiCentral)*std::sin(phi);
	double phi_flat = std::atan2(phi_numerator, phi_denominator);
	if(phi_flat < 0) phi_flat += M_PI*2.0;
	// double phi_flat_numerator = std::sin(theta)*std::sin(phi-m_phiCentral);
	// double phi_flat_denominator = std::cos(m_tilt)*std::sin(theta)*std::cos(phi-m_phiCentral) - std::sin(m_tilt)*std::cos(theta);
	// double phi_flat = std::atan2(phi_flat_numerator,phi_flat_denominator);
	// if(phi_flat < 0) phi_flat += 2*M_PI;

	//Calculate the *potential* R in the flat detector
	double r_numerator = m_translation.GetZ()*std::cos(phi)*std::sin(theta);
	double r_denominator = std::cos(phi_flat)*std::cos(m_phiCentral)*std::cos(m_tilt)*std::cos(theta) - std::sin(phi_flat)*std::sin(m_phiCentral)*std::cos(theta) - std::cos(phi_flat)*std::sin(m_tilt)*std::cos(phi)*std::sin(theta);
	double r_flat = r_numerator/r_denominator;

	//Calculate the distance from the origin to the hit on the detector
	//double R_to_detector = (r_flat*std::cos(phi_flat)*std::sin(m_tilt) + m_translation.GetZ())/std::cos(theta);

	//std::cout << "Theta = " << theta*rad2deg << ", phi = " << phi*rad2deg << std::endl;
	//std::cout << "r = " << r_flat << ", hitphi = " << phi_flat << std::endl << std::endl;


	//Check to see if our flat coords fall inside the flat detector
	if(IsInside(r_flat, phi_flat)) {
		int ringch, wedgech;
		if(phi_flat > M_PI) phi_flat -= 2.0*M_PI; //Need phi in terms of [-deltaPhi_flat/2, deltaPhi_flat/2]
		for(int i=0; i<m_nRings; i++) {
			if(IsRingTopEdge(r_flat, i) || IsRingBottomEdge(r_flat, i)) { //If it falls in the interstrip spacing, kill it
				ringch = -1;
				break;
			} else if(IsRing(r_flat, i)) {
				ringch = i;
				break;
			}
		}
		for(int i=0; i<m_nWedges; i++) {
			if(IsWedgeTopEdge(phi_flat, i) || IsWedgeBottomEdge(phi_flat, i)){ //If it falls in the interstrip spacing, kill it
				wedgech = -1;
				break;
			} else if(IsWedge(phi_flat, i)) {
				wedgech = i;
				break;
			}
		}
		return std::make_pair(ringch, wedgech);
	} else {
		return std::make_pair(-1,-1);
	}
}

/*
#include <cmath>
#include <utility>
#include "Vec3.h"
#include "Rotation.h"

std::pair<double, double> GetFlatCoordinatesFromOffsetRay(
    const Vec3& origin,
    double theta,
    double phi,
    double phiCentral,
    double tilt,
    double detectorZ
) {
    // Step 1: Define ray direction vector from spherical angles
    Vec3 direction;
    direction.SetVectorSpherical(1.0, theta, phi); // Unit vector

    // Step 2: Compute detector normal vector
    Vec3 detectorNormal(
        std::sin(tilt) * std::cos(phiCentral),
        std::sin(tilt) * std::sin(phiCentral),
        std::cos(tilt)
    );

    // Step 3: Compute detector position vector
    Vec3 detectorPos(0.0, 0.0, detectorZ);

    // Step 4: Compute intersection time t
    Vec3 delta = origin - detectorPos;
    double numerator = -delta.Dot(detectorNormal);
    double denominator = direction.Dot(detectorNormal);

    if (std::abs(denominator) < 1e-8) {
        // Ray is parallel to the detector plane
        return std::make_pair(-1.0, -1.0);
    }

    double t_hit = numerator / denominator;

    if (t_hit < 0) {
        // Hit is behind the ray origin
        return std::make_pair(-1.0, -1.0);
    }

    // Step 5: Compute hit point
    Vec3 hitPoint = Vec3(
        origin.GetX() + t_hit * direction.GetX(),
        origin.GetY() + t_hit * direction.GetY(),
        origin.GetZ() + t_hit * direction.GetZ()
    );

    // Step 6: Convert to detector-local coordinates
    Vec3 relativeHit = hitPoint - detectorPos;

    ZRotation rotZ(-phiCentral); // undo φ rotation
    YRotation rotY(-tilt);       // undo tilt

    // Combined rotation
    double temp[3];
    for (int i = 0; i < 3; ++i) {
        temp[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
            temp[i] += rotZ.m_matrix[i][j] * relativeHit.GetIndex(j);
        }
    }

    Vec3 intermediate(temp[0], temp[1], temp[2]);

    double local[3];
    for (int i = 0; i < 3; ++i) {
        local[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
            local[i] += rotY.m_matrix[i][j] * intermediate.GetIndex(j);
        }
    }

    Vec3 localHit(local[0], local[1], local[2]);

    // Step 7: Compute r_flat and phi_flat in local (detector) frame
    double x = localHit.GetX();
    double y = localHit.GetY();

    double r_flat = std::sqrt(x*x + y*y);
    double phi_flat = std::atan2(y, x);
    if (phi_flat < 0) phi_flat += 2.0 * M_PI;

    return std::make_pair(r_flat, phi_flat);
}
*/

/*
	Given a unit vector (R=1, theta, phi) which corresponds to some particle's trajectory
	and a vector containing translational offset information,
	determine whether that particle will intersect with this SABRE detector. If it does,
	determine which ring and wedge the hit occurs in.
*/
// std::pair<int, int> SABRE_Detector::GetOffsetTrajectoryRingWedge(double theta, double phi, const Vec3& offset){
// 	Vec3 direction(std::sin(theta)*std::cos(phi),
// 								 std::sin(theta)*std::sin(phi),
// 								 std::cos(theta));

// 	Vec3 planeNormal = GetNormTilted();
// 	Vec3 planePoint = m_translation;// GetHitCoordinates(8,4);//roughly somewhere near the center of the detector

// 	double denom = planeNormal.Dot(direction);
// 	if(std::fabs(denom) < 1e-6) return std::make_pair(-1,-1);

// 	double t = (planePoint - offset).Dot(planeNormal) / denom;
// 	if(t<0) return std::make_pair(-1,-1);

// 	Vec3 hitpoint = offset + direction*t;

// 	Vec3 shifted = hitpoint - m_translation;
// 	Vec3 detFrameVec = m_YRot.GetInverse()*(m_ZRot.GetInverse()*shifted);

// 	double r =std::sqrt(detFrameVec.GetX()*detFrameVec.GetX() + detFrameVec.GetY()*detFrameVec.GetY());
// 	double hitphi = std::atan2(detFrameVec.GetY(), detFrameVec.GetX());
// 	if(hitphi < 0) hitphi += 2.*M_PI;

// 	if(!IsInside(r, hitphi)) return std::make_pair(-1,-1);

// 	int ring = -1;
// 	int wedge = -1;

// 	if(hitphi > M_PI) hitphi -= 2.*M_PI;

// 	std::cout << "Theta = " << theta*rad2deg << ", phi = " << phi << std::endl;
// 	std::cout << "r = " << r << ", hitphi = " << hitphi << std::endl << std::endl;

// 	for(int i=0; i < m_nRings; i++){
// 		if(IsRingTopEdge(r,i) || IsRingBottomEdge(r,i)){
// 			return std::make_pair(-1,-1);
// 		}
// 		if(IsRing(r,i)){
// 			ring = i;
// 			break;
// 		}
// 	}

// 	for(int i=0; i<m_nWedges; i++){
// 		if(IsWedgeTopEdge(hitphi,i) || IsWedgeBottomEdge(hitphi,i)){
// 			return std::make_pair(-1,-1);
// 		}
// 		if(IsWedge(hitphi,i)){
// 			wedge = i;
// 			break;
// 		}
// 	}

// 	return std::make_pair(ring, wedge);

// }

std::pair<int,int> SABRE_Detector::GetOffsetTrajectoryRingWedge(double theta, double phi, const Vec3& offset){
	
	//establish direction vector of particle
	Vec3 direction;
	direction.SetVectorSpherical(1., theta, phi);

	Vec3 detectorNormal = GetNormTilted();
	Vec3 detectorPoint = m_translation;

	//parameterize intersection and calculate t
	Vec3 delta = offset - detectorPoint;
	double numerator = -delta.Dot(detectorNormal);
	double denominator = direction.Dot(detectorNormal);

	if(std::fabs(denominator) < 1e-6){
		return std::make_pair(-1,-1);//parallel
	}

	double t_hit = numerator/denominator;
	if(t_hit < 0){
		//hit is behind ray origin point 
		return std::make_pair(-1,-1);
	}

	Vec3 hitpoint = Vec3(offset.GetX() + t_hit*direction.GetX(),
						 offset.GetY() + t_hit*direction.GetY(),
						 offset.GetZ() + t_hit*direction.GetZ());
	Vec3 shifted = hitpoint - detectorPoint;
	Vec3 relativeHit = m_YRot.GetInverse()*(m_ZRot.GetInverse()*shifted);

	double x = relativeHit.GetX();
	double y = relativeHit.GetY();

	double r_flat = std::sqrt(x*x + y*y);
	double phi_flat = std::atan2(y,x);
	if(phi_flat < 0) phi_flat += M_PI*2.;

	if(IsInside(r_flat, phi_flat)){
		int ringch, wedgech;
		if(phi_flat > M_PI) phi_flat -= 2.*M_PI;
		for(int i=0; i<m_nRings; i++){
			if(IsRingTopEdge(r_flat,i) || IsRingBottomEdge(r_flat,i)){
				ringch = -1;
				break;
			} else if(IsRing(r_flat,i)){
				ringch = i;
				break;
			}
		}
		for(int i=0; i<m_nWedges; i++){
			if(IsWedgeTopEdge(phi_flat,i) || IsWedgeBottomEdge(phi_flat,i)){
				wedgech = -1;
				break;
			} else if(IsWedge(phi_flat,i)){
				wedgech = i;
				break;
			}
		}
		return std::make_pair(ringch,wedgech);
	} else {
		return std::make_pair(-1,-1);
	}
}


/*
	Given a ring/wedge of this SABRE detector, calculate the coordinates of a hit.
	Currently gives a point in the *center* of the pixel. Better method would be to
	randomly wiggle the point within the pixel. Method intended for use with data, or
	to smear out simulated data to mimic real data.
*/
Vec3 SABRE_Detector::GetHitCoordinates(int ringch, int wedgech) {
	if(!CheckRingChannel(ringch) || !CheckWedgeChannel(wedgech)) {
		return Vec3();
	}

	double r_center  = m_Rinner + (0.5+ringch)*m_deltaR_flat_ring;
	double phi_center = -m_deltaPhi_flat/2.0 + (0.5+wedgech)*m_deltaPhi_flat_wedge;
	double x = r_center*std::cos(phi_center);
	double y = r_center*std::sin(phi_center);
	double z = 0;

	Vec3 hit(x, y, z);

	return TransformToTiltedFrame(hit);
}
/*
	Given a ring/wedge of this SABRE detector, calculate the coordinates of a hit.
	Adds a "random wiggle" of the point in the area of the pixel. This avoids
	dense colletion of points near the edges when just randomly sampling r, phi separately.
*/
//added by zmp:
Vec3 SABRE_Detector::GetHitCoordinatesRandomWiggle(int ringch, int wedgech){
	if(!CheckRingChannel(ringch) || !CheckWedgeChannel(wedgech)){
		return Vec3();
	}

	//define pixel boundaries
	double r_min = m_Rinner + ringch*m_deltaR_flat_ring;
	double r_max = m_Rinner + (ringch + 1)*m_deltaR_flat_ring;
	double phi_min = -m_deltaPhi_flat/2.0 + wedgech*m_deltaPhi_flat_wedge;
	double phi_max = -m_deltaPhi_flat/2.0 + (wedgech + 1)*m_deltaPhi_flat_wedge;

	//random number generation
	static thread_local std::mt19937 gen(std::random_device{}());
	std::uniform_real_distribution<double> ur(0.,1.);

	//uniform area sampling
	double r = std::sqrt(ur(gen)*(r_max*r_max - r_min*r_min) + r_min*r_min);
	double phi = phi_min + ur(gen)*(phi_max - phi_min);

	double x = r*std::cos(phi);
	double y = r*std::sin(phi);
	double z = 0.;

	Vec3 hit(x,y,z);

	return TransformToTiltedFrame(hit);
}

void SABRE_Detector::WriteTransformedCorners(std::ofstream& outfile) {
	//outfile << std::setprecision(10);

	for(int ring=0; ring<m_nRings; ring++){
		for(int wedge=0; wedge<m_nWedges; wedge++){
			double r_inner = m_Rinner + ring*m_deltaR_flat_ring;
			double r_outer = m_Rinner + (ring+1)*m_deltaR_flat_ring;
			double phi_start = -m_deltaPhi_flat/2. + wedge*m_deltaPhi_flat_wedge;
			double phi_end = -m_deltaPhi_flat/2. + (wedge+1)*m_deltaPhi_flat_wedge;

			Vec3 corner0_flat(r_outer*std::cos(phi_start),r_outer*std::sin(phi_start),0);
			Vec3 corner1_flat(r_inner*std::cos(phi_start),r_inner*std::sin(phi_start),0);
			Vec3 corner2_flat(r_inner*std::cos(phi_end),r_inner*std::sin(phi_end),0);
			Vec3 corner3_flat(r_outer*std::cos(phi_end),r_outer*std::sin(phi_end),0);

			Vec3 corner0_tilt = TransformToTiltedFrame(corner0_flat);
			Vec3 corner1_tilt = TransformToTiltedFrame(corner1_flat);
			Vec3 corner2_tilt = TransformToTiltedFrame(corner2_flat);
			Vec3 corner3_tilt = TransformToTiltedFrame(corner3_flat);

			outfile << ring << " " << wedge << " "
						  << corner0_tilt.GetX() << " " << corner0_tilt.GetY() << " "
						  << corner1_tilt.GetX() << " " << corner1_tilt.GetY() << " " 
						  << corner2_tilt.GetX() << " " << corner2_tilt.GetY() << " "
						  << corner3_tilt.GetX() << " " << corner3_tilt.GetY() << "\n";
						  //<< ring << " " << wedge << "\n";
		}
	}
}
using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>

#include "SABRE_Detector.h"
#include "Vec3.h"
#include "Rotation.h"
#include "Beamspot.h"
#include "SABRE_EnergyResolutionModel.h"

static const std::pair<int, int> offsets[] = {
	{112,40},	//detector0
	{96,32},	//detector1
	{80,16},	//detector2
	{64,24},	//detector3
	{48,0}		//detector4
};

int main(int argc, char * argv[]){
	int eoev = -1;
	int nevents = 0;
	bool detected1 = false;
	bool detected2 = false;
	//bool detectedboth = false;
	int hit1 = 0;
	int hit2 = 0;
	int hitboth = 0;
	float dx = 0, dy = 0;
	float e1, theta1, phi1, thetacm, e2, theta2, phi2;//kin2mc output

	//set seed for RNG
	srand(time(NULL));

	for(int narg=0; narg<argc; narg++){
		cout << "argument " << narg << " = " << argv[narg] << endl;
	}

	if(argc<3){
		cerr << "Error: Please provide input and output file name as command-line arguments!" << endl;
		return 1;
	}

	ifstream infile(argv[1]);
	ofstream outfile(argv[2]);

	cout << "Enter beamspot dx, dy (1/line):" << endl;
	cin >> dx;
	cin >> dy;

	Beamspot beamspot;
	beamspot.SetXMax(dx);
	beamspot.SetYMax(dy);
	beamspot.SetXOffset(0);
	beamspot.SetYOffset(0);

	cout << endl;
	cout << "Processing physics data file " << argv[1] << endl;
	cout << "Writing to output file " << argv[2] << endl << endl;
	cout << "Beam spot spread is (dx,dy) = (" << dx << ", " << dy << ")" << endl << endl;

	double DEG2RAD = M_PI/180.;
	double INNER_R = 0.0326;
	double OUTER_R = 0.1351;
	double TILT = 40.;
	double ZDIST = -0.1245;
	double PHI_COVERAGE = 54.4;
	vector<double> PHI = {306.,18.,234.,162.,90.};//in order of detectors 0-4 from Rachel sabre powerpoint
	vector<SABRE_Detector*> SABRE_Array;
	vector<int> SABRE_Array_hits = {0,0,0,0,0};
	for(size_t i=0; i<PHI.size(); i++){
		SABRE_Array.push_back(new SABRE_Detector(INNER_R, OUTER_R, PHI_COVERAGE*DEG2RAD, PHI[i]*DEG2RAD, TILT*DEG2RAD, ZDIST));
		//cout << "Successfully created SABRE_Detector at PHI[" << i << "] = " << PHI[i] << endl;
	}
	cout << endl;

	//prepare the SABRE_EnergyResolutionModels here (1/SABRE detector, so 5 total in the array)
	std::vector<SABRE_EnergyResolutionModel*> SABREARRAY_EnergyResolutionModels;
	for(size_t i=0; i<SABRE_Array.size(); i++){
		SABREARRAY_EnergyResolutionModels.push_back(new SABRE_EnergyResolutionModel(0.05,1.));////sigma of 0.05 MeV for all rings/wedges and threshold of 1MeV for all rings/wedges -- default, but can update with setters and eventually read in from a file!
	}

/*want to debug this so that I can actually use beamspot, but for now let's just assume perfect precision beamspot at (0,0,0)
	while(infile >> e1 >> theta1 >> phi1 >> thetacm >> e2 >> theta2 >> phi2){
		nevents += 1;
		detected1 = false;
		detected2 = false;
		//detectedboth = false;

		beamspot.Spread();

		outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << thetacm << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << endl;
		if(nevents%50000==0) cout << "Processed " << nevents << " events..." << endl;


		for(size_t i=0; i<SABRE_Array.size(); i++){

			//check the first ejectile:
			if(beamspot.Set(SABRE_Array[i]->GetRingTiltCoords(0,0), SABRE_Array[i]->GetNormTilted(), theta1, phi1)){
				//cout << "true" << endl;
				//this means that the trajectory of the particle intersects SABRE_Array[i] when eminating from an origin of (Dx, Dy, 0) instead of (0,0,0)
				float newTheta, newPhi;
				newTheta = beamspot.GetTheta();
				newPhi = beamspot.GetPhi();
				//since we know that the trajectory of the particle intersects SABRE_Array[i], we can just pass these new theta and new phi to the original function to get ring,wedge!
				pair<int,int> hit1_rw = SABRE_Array[i]->GetTrajectoryRingWedge(newTheta*DEG2RAD,newPhi*DEG2RAD);
				if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){
					//cout << "true" << endl;
					Vec3 localCoords = SABRE_Array[i]->GetHitCoordinates(hit1_rw.first,hit1_rw.second);
					double smeared_ering, smeared_ewedge;
					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit1_rw.first,e1,smeared_ering) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit1_rw.second,e1,smeared_ewedge)){
						outfile << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smeared_ering << "\t" << smeared_ewedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
						detected1 = true;
						hit1 += 1;
						SABRE_Array_hits[i] += 1;
					}
				}
			}

			//check the second ejectile:
			if(beamspot.Set(SABRE_Array[i]->GetRingTiltCoords(0,0), SABRE_Array[i]->GetNormTilted(), theta2, phi2)){
				//this means that the trajectory of the particle intersects SABRE_Array[i]
				float newTheta, newPhi;
				newTheta = beamspot.GetTheta();
				newPhi = beamspot.GetPhi();
				//since we know the trajectory of the particle intersects SABRE_Array[i], we can just pass these new theta and new phi to the original function to get ring,wedge!
				pair<int,int> hit2_rw = SABRE_Array[i]->GetTrajectoryRingWedge(newTheta*DEG2RAD,newPhi*DEG2RAD);
				if(hit2_rw.first != -1 && hit2_rw.second != -1 && !detected2){
					Vec3 localCoords = SABRE_Array[i]->GetHitCoordinates(hit2_rw.first,hit2_rw.second);
					double smeared_ering, smeared_ewedge;
					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit2_rw.first,e2,smeared_ering) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit2_rw.second,e2,smeared_ewedge)){
						outfile << 200+i << "\t" << hit2_rw.first << "\t" << hit2_rw.second << "\t" << smeared_ering << "\t" << smeared_ewedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
						detected2 = true;
						hit2 += 1;
						SABRE_Array_hits[i] += 1;
					}
				}
			}

			if(detected1 && detected2){
				hitboth += 1;
				break;
			}

		}

		outfile << eoev << endl;
	}
*/

	while(infile >> e1 >> theta1 >> phi1 >> thetacm >> e2 >> theta2 >> phi2){
		nevents += 1;
		detected1 = false;
		detected2 = false;

		//beamspot.Spread();//not necessary right now since Set() not working right

		outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << thetacm << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << endl;
		if(nevents%50000==0) cout << "Processed " << nevents << " events..." << endl;

		for(size_t i=0; i<SABRE_Array.size(); i++){
			double smearedERing, smearedEWedge;

			pair<int,int> hit1_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta1*DEG2RAD, phi1*DEG2RAD);
			if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){
				if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit1_rw.first,e1,smearedERing), SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit1_rw.second,e1,smearedEWedge)){
					Vec3 localCoords = SABRE_Array[i]->GetHitCoordinates(hit1_rw.first,hit1_rw.second);
					outfile << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
					detected1 = true;
					hit1+=1;
					SABRE_Array_hits[i] += 1;
				}
			}

			smearedERing = 0.;
			smearedEWedge = 0.;
			pair<int,int> hit2_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta2*DEG2RAD,phi2*DEG2RAD);
			if(hit2_rw.first != -1 && hit2_rw.second != -1 && !detected2){
				if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit2_rw.first,e2,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit2_rw.second,e2,smearedEWedge)){
					Vec3 localCoords = SABRE_Array[i]->GetHitCoordinates(hit2_rw.first,hit2_rw.second);
					outfile << 200+i << "\t" << hit2_rw.first << "\t" << hit2_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
					detected2 = true;
					hit2+=1;
					SABRE_Array_hits[i] += 1;
				}
			}
		}

		outfile << eoev << endl;

	}

	infile.close();
	outfile.close();

	cout << endl;
	cout << "Files are closed. Processed " << nevents << " events." << endl;
	cout << "Events with ejectile in SABRE: " << hit1 << endl;
	cout << "Events with recoil in SABRE: " << hit2 << endl;
	cout << "Events with both ejectile and recoil in SABRE: " << hitboth << endl << endl;
	for(int i=0; i<5; i++){
		cout << "Detector_" << i << " had total hits = " << SABRE_Array_hits[i] << endl;
	}
	cout << endl;
	cout << "Output can be found here: " << argv[2] << endl;
	cout << endl;

	for(SABRE_Detector* detector : SABRE_Array){
		delete detector;
	}
	SABRE_Array_hits.clear();
}

/*Geometric Efficiencies Commented Out Below*/
/*
	double DEG2RAD = M_PI/180.;
	double INNER_R = 0.0326;
	double OUTER_R = 0.1351;
	double TILT = 40.;
	double ZDIST = -0.1245;
	double PHI_COVERAGE = 54.4;
	vector<double> PHI = {234.,162.,306.,18.,90.};
	//time to test if we can make a single sabre detector:
	vector<SABRE_Detector*> SABRE_Array;
	vector<int> SABRE_Array_hits = {0,0,0,0,0};

	vector<double> SABRE_Array_minTheta = {181.,181.,181.,181.,181.};
	vector<double> SABRE_Array_maxTheta = {-1.,-1.,-1.,-1.,-1.};
	vector<double> SABRE_Array_minPhi = {361.,361.,361.,361.,361.};
	vector<double> SABRE_Array_maxPhi = {-1.,-1.,-1.,-1.,-1.};
	//SABRE_Array.push_back(new SABRE_Detector(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI[0]*DEG2RAD,TILT*DEG2RAD,ZDIST));
	for(size_t i=0; i<PHI.size(); i++){
		SABRE_Array.push_back(new SABRE_Detector(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI[i]*DEG2RAD,TILT*DEG2RAD,ZDIST));
		cout << "Sucessfully created SABRE_Detector at PHI[" << i << "] = " << PHI[i] << endl;
	}

	double thetastep = 0.1, phistep = 0.1;
	int totalrays = 0, hits = 0;
	cout << "Casting " << (180./thetastep)*(360./phistep) << " rays to determine geometric efficencies..." << endl;
	for(double theta=0.; theta<=180.; theta+=thetastep){
		for(double phi=0.; phi<360.; phi+=phistep){
			totalrays += 1;
			for(int i=0; i<5; i++){
				Vec3 intersection = SABRE_Array[i]->GetTrajectoryCoordinates(theta*DEG2RAD,phi*DEG2RAD);
				if(intersection.GetX() != 0 && intersection.GetY() != 0 && intersection.GetZ() != 0){
					hits += 1;
					SABRE_Array_hits[i] += 1;
					if(SABRE_Array_minTheta[i] > theta) SABRE_Array_minTheta[i] = theta;
					if(SABRE_Array_maxTheta[i] < theta) SABRE_Array_maxTheta[i] = theta;
					if(SABRE_Array_minPhi[i] > phi) SABRE_Array_minPhi[i] = phi;
					if(SABRE_Array_maxPhi[i] < phi) SABRE_Array_maxPhi[i] = phi;
					break;
				}
			}
		}
	}

	cout << "Total rays cast: " << totalrays << endl;
	cout << "Total rays hit: " << hits << endl;
	cout << "Geometric Efficiency: " << 100.*static_cast<double>(hits)/static_cast<double>(totalrays) << endl << endl;
	for(size_t i=0; i<SABRE_Array.size(); i++){
		cout << "Detector " << i << " hits = " << SABRE_Array_hits[i] << endl;
		cout << "\tMin theta = " << SABRE_Array_minTheta[i] << endl;
		cout << "\tMax theta = " << SABRE_Array_maxTheta[i] << endl << endl;
		cout << "\tMin Phi = " << SABRE_Array_minPhi[i] << endl;
		cout << "\tMax Phi = " << SABRE_Array_maxPhi[i] << endl << endl;
	}
*/



//}
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
	bool detected3 = false;
	bool detected4 = false;
	//bool detectedboth = false;
	int hit1 = 0;
	int hit2 = 0;
	int hit3 = 0;
	int hit4 = 0;
	int hitboth = 0;
	//float dx = 0, dy = 0;
	float e1, theta1, phi1, thetacm, e2, theta2, phi2;//kin2mc output
	float e3, theta3, phi3, e4, theta4, phi4;//kin3/4mc output

	//set seed for RNG
	srand(time(NULL));

	// for(int narg=0; narg<argc; narg++){
	// 	cout << "argument " << narg << " = " << argv[narg] << endl;
	// }

	if(argc!=4){
		cerr << "Error: Please provide input as command line arguments!" << endl;
		cerr << "Expected arguments in order of: X kinInputFile.out detOutputFile.det" << endl;
		cerr << "Where X = 2,3,4 for kin2mc, kin3mc, kin4mc input files" << endl;
		cerr << "And kinInputFile.out is the filename (or path to) the kinXmc output file" << endl;
		cerr << "And detOutputFile.det is the filename (or path to) the detection output file (what this code writes to)" << endl;
		cout << endl;
		return 1;
	}

	int kinX = std::stoi(argv[1]);
	if(kinX != 2 && kinX != 3 && kinX != 4){
		cerr << "Error: Invalid kinematics type. Use 2 for 2-body, 3 for 3-body, 4 for 4-body!" << endl;
		return 1;
	}

	ifstream infile(argv[2]);
	ofstream outfile(argv[3]);

	// cout << "Enter beamspot dx, dy (1/line):" << endl;
	// cin >> dx;
	// cin >> dy;

	// Beamspot beamspot;
	// beamspot.SetXMax(dx);
	// beamspot.SetYMax(dy);
	// beamspot.SetXOffset(0);
	// beamspot.SetYOffset(0);

	cout << endl;
	cout << "Kin" << kinX << "mc selected!" << endl;
	cout << "Processing physics data file " << argv[2] << endl;
	cout << "Writing to output file " << argv[3] << endl << endl;
	//cout << "Beam spot spread is (dx,dy) = (" << dx << ", " << dy << ")" << endl << endl;

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

	int onePartHits = 0, twoPartHits=0, threePartHits=0, fourPartHits=0;

	//need the transformed corners for each detector, so lets just output that here really quickly:
	// ofstream cornerfiles[5];
	// cornerfiles[0].open("corners/SABRE0_corners.txt");
	// cornerfiles[1].open("corners/SABRE1_corners.txt");
	// cornerfiles[2].open("corners/SABRE2_corners.txt");
	// cornerfiles[3].open("corners/SABRE3_corners.txt");
	// cornerfiles[4].open("corners/SABRE4_corners.txt");

	// for(size_t i=0; i<SABRE_Array.size(); i++){
	// 	SABRE_Detector* detector = SABRE_Array[i];
	// 	//cornerfiles[i] << "SABRE_" << i << endl;
	// 	detector->WriteTransformedCorners(cornerfiles[i]);
	// 	cornerfiles[i] << endl;
		
	// }

	// for(int i=0; i<5; i++){
	// 	cornerfiles[i].close();
	// }

	//prepare the SABRE_EnergyResolutionModels here (1/SABRE detector, so 5 total in the array)
	std::vector<SABRE_EnergyResolutionModel*> SABREARRAY_EnergyResolutionModels;
	for(size_t i=0; i<SABRE_Array.size(); i++){
		SABREARRAY_EnergyResolutionModels.push_back(new SABRE_EnergyResolutionModel(0.05,1.));////sigma of 0.05 MeV for all rings/wedges and threshold of 1MeV for all rings/wedges -- default, but can update with setters and eventually read in from a file!
	}

	if(kinX == 2){//kin2mc
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
					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit1_rw.first,e1,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit1_rw.second,e1,smearedEWedge)){
						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
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
						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit2_rw.first,hit2_rw.second);
						outfile << 200+i << "\t" << hit2_rw.first << "\t" << hit2_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
						detected2 = true;
						hit2+=1;
						SABRE_Array_hits[i] += 1;
					}
				}
			}

			outfile << eoev << endl;

		}
	} else if(kinX == 3){//kin3mc
		while(infile >> e1 >> theta1 >> phi1 >> e2 >> theta2 >> phi2 >> e3 >> theta3 >> phi3 >> e4 >> theta4 >> phi4){
			nevents += 1;
			detected1 = false;
			detected3 = false;
			detected4 = false;

			//outfile << std::format("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",e1,theta1,phi1,e2,theta2,phi2,e3,theta3,phi3,e4,theta4,phi4) << endl;
			outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\t" << e3 << "\t" << theta3 << "\t" << phi3 << "\t" << e4 << "\t" << theta4 << "\t" << phi4 << endl;
			if(nevents%50000==0) cout << "Processed " << nevents << " events..." << endl;

			for(size_t i=0; i<SABRE_Array.size(); i++){
				double smearedERing, smearedEWedge;

				pair<int,int> hit1_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD);
				if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){
					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit1_rw.first,e1,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit1_rw.second,e1,smearedEWedge)){
						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
						//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",100+i,hit1_rw.first,hit1_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
						outfile << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
						detected1=true;
						hit1+=1;
						SABRE_Array_hits[i] += 1;
					}
				}

				smearedERing = 0.;
				smearedEWedge = 0.;
				pair<int,int> hit3_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta3*DEG2RAD,phi3*DEG2RAD);
				if(hit3_rw.first != -1 && hit3_rw.second != -1 && !detected2){
					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit3_rw.first,e3,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit3_rw.second,e3,smearedEWedge)){
						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit3_rw.first,hit3_rw.second);
						//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",300+i,hit3_rw.first,hit3_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
						outfile << 300+i << "\t" << hit3_rw.first << "\t" << hit3_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
						detected3 = true;
						hit3 += 1;
						SABRE_Array_hits[i] += 1;
					}
				}

				smearedERing = 0.;
				smearedEWedge = 0.;
				pair<int,int> hit4_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta4*DEG2RAD,phi4*DEG2RAD);
				if(hit4_rw.first != -1 && hit4_rw.second != -1 && !detected4){
					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit4_rw.first,e4,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit4_rw.second,e4,smearedEWedge)){
						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit4_rw.first,hit4_rw.second);
						//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",400+i,hit4_rw.first,hit4_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
						outfile << 400+i << "\t" << hit4_rw.first << "\t" << hit4_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
						detected4 = true;
						hit4 += 1;
						SABRE_Array_hits[i] += 1;
					}
				}
			}

			outfile << eoev << endl;

			if((detected1&&!detected3&&!detected4) || (!detected1&&detected3&&!detected4) || (!detected1&&!detected3&&detected4)){
				onePartHits += 1;
			} else if((detected1&&detected3&&!detected4) || (!detected1&&detected3&&detected4) || (detected1&&!detected3&&detected4)){
				twoPartHits += 1;
			} else if((detected1&&detected3&&detected4)){
				threePartHits += 1;
			}
		}
	} else if(kinX == 4){//kin4mc
		while(infile >> e1 >> theta1 >> phi1 >> e2 >> theta2 >> phi2 >> e3 >> theta3 >> phi3 >> e4 >> theta4 >> phi4){
			nevents += 1;
			detected1 = false;
			detected2 = false;
			detected3 = false;
			detected4 = false;

			outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\t" << e3 << "\t" << theta3 << "\t" << phi3 << "\t" << e4 << "\t" << theta4 << "\t" << phi4 << endl;
			if(nevents%50000==0) cout << "Processed " << nevents << " events..." << endl;

			for(size_t i=0; i<SABRE_Array.size(); i++){
				double smearedERing, smearedEWedge;

				pair<int,int> hit1_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD);
				if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){
					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit1_rw.first,e1,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit1_rw.second,e1,smearedEWedge)){
						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
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
						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit2_rw.first,hit2_rw.second);	
						outfile << 200+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
						detected2 = true;
						hit2+=1;
						SABRE_Array_hits[i] += 1;
					}
				}

				smearedERing = 0.;
				smearedEWedge = 0.;
				pair<int,int> hit3_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta3*DEG2RAD,phi3*DEG2RAD);
				if(hit3_rw.first != -1 && hit3_rw.second != -1 && !detected3){
					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit3_rw.first,e3,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit3_rw.second,e3,smearedEWedge)){
						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit3_rw.first,hit3_rw.second);
						outfile << 300+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
						detected3 = true;
						hit3+=1;
						SABRE_Array_hits[i] += 1;
					}
				}

				smearedERing = 0.;
				smearedEWedge = 0.;
				pair<int,int> hit4_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta4*DEG2RAD,phi4*DEG2RAD);
				if(hit4_rw.first != -1 && hit4_rw.second != -1 && !detected4){
					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit4_rw.first,e4,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit4_rw.second,e4,smearedEWedge)){
						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit4_rw.first,hit4_rw.second);
						outfile << 400+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
						detected4 = true;
						hit4+=1;
						SABRE_Array_hits[i] += 1;
					}
				}
			}

			outfile << eoev << endl;

			if((detected1&&!detected2&&!detected3&&!detected4) || (!detected1&&detected2&&!detected3&&!detected4) || (!detected1&&!detected2&&detected3&&!detected4) || (!detected1&&!detected2&&!detected3&&detected4)){
				onePartHits += 1;
			} else if((detected1&&detected2&&!detected3&&!detected4) || (!detected1&&detected2&&detected3&&!detected4) || (!detected1&&!detected2&&detected3&&detected4) || (detected1&&!detected2&&!detected3&&detected4)){
				twoPartHits += 1;
			} else if((detected1&&detected2&&detected3&&!detected4) || (!detected1&&detected2&&detected3&&detected4) || (detected1&&!detected2&&detected3&&detected4) || (detected1&&detected2&&!detected3&&detected4)){
				threePartHits += 1;
			} else if((detected1&&detected2&&detected3&&detected4)){
				fourPartHits += 1;
			}
		}
	}

	infile.close();
	outfile.close();

	cout << endl;
	cout << "Files are closed. Processed " << nevents << " events." << endl;
	if(kinX==2||kinX==3||kinX==4) cout << "Events with ejectile in SABRE: " << hit1 << endl;
	if(kinX==2) cout << "Events with recoil in SABRE: " << hit2 << endl;
	if(kinX==3) cout << "Events with bu1 in SABRE: " << hit3 << "\nEvents with bu2 in SABRE: " << hit4 << endl;
	if(kinX==4) cout << "Events with bu1 in SABRE: " << hit2 << "\nEvents with bu2 in SABRE: " << hit3 << "\nEvents with bu3 in SABRE: " << hit4 << endl; 
	//cout << "Events with both ejectile and recoil in SABRE: " << hitboth << endl << endl;
	if(onePartHits>0) cout << "1 Particle Events: " << onePartHits << endl;
	if(twoPartHits>0) cout << "2 Particle Events: " << twoPartHits << endl;
	if(threePartHits>0) cout << "3 Particle Events: " << threePartHits << endl;
	if(fourPartHits>0) cout << "4 Particle Events: " << fourPartHits << endl;
	for(int i=0; i<5; i++){
		cout << "Detector_" << i << " had total hits = " << SABRE_Array_hits[i] << endl;
	}
	cout << endl;
	cout << "Output can be found here: " << argv[3] << endl;
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
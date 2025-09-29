using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include "ConsoleColorizer.h"
#include "SABREsim.h"

static const std::pair<int, int> offsets[] = {
	{112,40},	//detector0
	{96,32},	//detector1
	{80,16},	//detector2
	{64,24},	//detector3
	{48,0}		//detector4
};

//old main - keep for now, but once confirmed working will remove and re-commit
// int main(int argc, char * argv[]){
// 	int eoev = -1;
// 	int nevents = 0;
// 	bool detected1 = false;
// 	bool detected2 = false;
// 	bool detected3 = false;
// 	bool detected4 = false;
// 	//bool detectedboth = false;
// 	int hit1 = 0;//num events where particle 1 is seen in SABRE
// 	int hit2 = 0;//num events where particle 2 is seen in SABRE
// 	int hit3 = 0;//num events where particle 3 is seen in SABRE
// 	int hit4 = 0;//num events where particle 4 is seen in SABRE
// 	int hit34 = 0;//num events where particles 3 and 4 are seen in SABRE in coincidence
// 	int hitonly3 = 0;//num events where only particle 3 is detected (no particle 4 seen)
// 	int hitonly4 = 0;//num events where only particle 4 is detected (no particle 3 seen)
// 	int hit1only = 0;//num events in kin2mc where only part1 is detected in SABRE (ej)
// 	int hit2only = 0;//num events in kin2mc where only part2 is detected in SABRE (rec)
// 	int hitboth = 0;//num events in kin2mc where both ej and rec detected in SABRE
// 	//int hitboth = 0;
// 	//float dx = 0, dy = 0;
// 	float e1, theta1, phi1, thetacm, e2, theta2, phi2;//kin2mc output
// 	float e3, theta3, phi3, e4, theta4, phi4;//kin3/4mc output

// 	//set seed for RNG
// 	srand(time(NULL));

// 	// for(int narg=0; narg<argc; narg++){
// 	// 	cout << "argument " << narg << " = " << argv[narg] << endl;
// 	// }

// 	if(argc == 2 && (std::string(argv[1]) == "help" || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")){
// 		ConsoleColorizer::PrintWhite("Usage:\n\t");
// 		ConsoleColorizer::PrintYellow("./SABREsim X kinXmcFile.out detOutputFile.det\n\n");
// 		ConsoleColorizer::PrintWhite("Where ");
// 		ConsoleColorizer::PrintYellow("X = 2,3,4 ");
// 		ConsoleColorizer::PrintWhite("for kin2mc, kin3mc, kin4mc input  files\nAnd ");
// 		ConsoleColorizer::PrintYellow("kinInputFile.out ");
// 		ConsoleColorizer::PrintWhite("is the filename (or path to) the kinXmc output file\nAnd ");
// 		ConsoleColorizer::PrintYellow("detOutputFile.det ");
// 		ConsoleColorizer::PrintWhite("is the filename (or path to) the detection output file (what this code writes to)\n\nExample command:\n");
// 		ConsoleColorizer::PrintYellow("\t./SABREsim 2 ../kinmc/TEST.out ../det/TEST.det\n\n");
// 		return 1;
// 	}

// 	if(argc!=4){
// 		// cerr << "Error: Please provide input as command line arguments!" << endl;
// 		// cerr << "Expected arguments in order of: X kinInputFile.out detOutputFile.det" << endl;
// 		// cerr << "Where X = 2,3,4 for kin2mc, kin3mc, kin4mc input files" << endl;
// 		// cerr << "And kinInputFile.out is the filename (or path to) the kinXmc output file" << endl;
// 		// cerr << "And detOutputFile.det is the filename (or path to) the detection output file (what this code writes to)" << endl;
// 		// cout << endl;
// 		ConsoleColorizer::PrintRed("\nError! Please provide input as command line arguments!\nExpected arguments in order of: ");
// 		ConsoleColorizer::PrintYellow("X kinXmcFile.out detOutputFile.det\n");
// 		ConsoleColorizer::PrintRed("Where ");
// 		ConsoleColorizer::PrintYellow("X = 2,3,4 ");
// 		ConsoleColorizer::PrintRed("for kin2mc, kin3mc, kin4mc input  files\nAnd ");
// 		ConsoleColorizer::PrintYellow("kinInputFile.out ");
// 		ConsoleColorizer::PrintRed("is the filename (or path to) the kinXmc output file\nAnd ");
// 		ConsoleColorizer::PrintYellow("detOutputFile.det ");
// 		ConsoleColorizer::PrintRed("is the filename (or path to) the detection output file (what this code writes to)\n\nExample command:\n");
// 		ConsoleColorizer::PrintYellow("\t./SABREsim 2 ../kinmc/TEST.out ../det/TEST.det\n\n");
// 		return 1;
// 	}

// 	int kinX = std::stoi(argv[1]);
// 	if(kinX != 2 && kinX != 3 && kinX != 4){
// 		cerr << "Error: Invalid kinematics type. Use 2 for 2-body, 3 for 3-body, 4 for 4-body!" << endl;
// 		return 1;
// 	}

// 	ifstream infile(argv[2]);
// 	ofstream outfile(argv[3]);

// 	// cout << "Enter beamspot dx, dy (1/line):" << endl;
// 	// cin >> dx;
// 	// cin >> dy;

// 	// Beamspot beamspot;
// 	// beamspot.SetXMax(dx);
// 	// beamspot.SetYMax(dy);
// 	// beamspot.SetXOffset(0);
// 	// beamspot.SetYOffset(0);

// 	cout << endl;
// 	cout << "Kin" << kinX << "mc selected!" << endl;
// 	cout << "Processing physics data file " << argv[2] << endl;
// 	cout << "Writing to output file " << argv[3] << endl << endl;
// 	//cout << "Beam spot spread is (dx,dy) = (" << dx << ", " << dy << ")" << endl << endl;

// 	double DEG2RAD = M_PI/180.;
// 	double INNER_R = 0.0326;
// 	double OUTER_R = 0.1351;
// 	double TILT = 40.;
// 	double ZDIST = -0.1245;
// 	double PHI_COVERAGE = 54.4;
// 	vector<double> PHI = {306.,18.,234.,162.,90.};//in order of detectors 0-4 from Rachel sabre powerpoint
// 	vector<SABRE_Detector*> SABRE_Array;
// 	vector<int> SABRE_Array_hits = {0,0,0,0,0};
// 	for(size_t i=0; i<PHI.size(); i++){
// 		SABRE_Array.push_back(new SABRE_Detector(INNER_R, OUTER_R, PHI_COVERAGE*DEG2RAD, PHI[i]*DEG2RAD, TILT*DEG2RAD, ZDIST));
// 		cout << "Successfully created SABRE_Detector at PHI[" << i << "] = " << PHI[i] << endl;
// 		cout << "\tNormal X: " << SABRE_Array[i]->GetNormTilted().GetX() << endl;
// 		cout << "\tNormal Y: " << SABRE_Array[i]->GetNormTilted().GetY() << endl;
// 		cout << "\tNormal Z: " << SABRE_Array[i]->GetNormTilted().GetZ() << endl;
// 		cout << endl;
// 	}
// 	cout << endl;

// 	int onePartHits = 0, twoPartHits=0, threePartHits=0, fourPartHits=0;

// 	//need the transformed corners for each detector, so lets just output that here really quickly:
// 	// ofstream cornerfiles[5];
// 	// cornerfiles[0].open("corners/SABRE0_corners.txt");
// 	// cornerfiles[1].open("corners/SABRE1_corners.txt");
// 	// cornerfiles[2].open("corners/SABRE2_corners.txt");
// 	// cornerfiles[3].open("corners/SABRE3_corners.txt");
// 	// cornerfiles[4].open("corners/SABRE4_corners.txt");

// 	// for(size_t i=0; i<SABRE_Array.size(); i++){
// 	// 	SABRE_Detector* detector = SABRE_Array[i];
// 	// 	//cornerfiles[i] << "SABRE_" << i << endl;
// 	// 	detector->WriteTransformedCorners(cornerfiles[i]);
// 	// 	cornerfiles[i] << endl;
		
// 	// }

// 	// for(int i=0; i<5; i++){
// 	// 	cornerfiles[i].close();
// 	// }

// 	//prepare the SABRE_EnergyResolutionModels here (1/SABRE detector, so 5 total in the array)
// 	//REWRITE this to read resolution in from a file, much easier than re-making everytime for a change in sigma!
// 	std::vector<SABRE_EnergyResolutionModel*> SABREARRAY_EnergyResolutionModels;
// 	for(size_t i=0; i<SABRE_Array.size(); i++){
// 		SABREARRAY_EnergyResolutionModels.push_back(new SABRE_EnergyResolutionModel(0.050,0.1));////sigma of 0.050 MeV for all rings/wedges and threshold of 0.100 MeV for all rings/wedges -- default, but can update with setters and eventually read in from a file!
// 	}

// 	//set up TargetEnergyLosses here:
// 	//6Li in LiF target (75 ug/cm^2) for E(6Li) from 0 keV to 20000 keV
// 	std::string targetEnergyLossPath_6Li_in_LiF = "../config/TargetELoss_6Li_in_LiF.conf";
// 	TargetEnergyLoss* targetLoss_6Li_in_LiF = TargetEnergyLoss::LoadFromConfigFile(targetEnergyLossPath_6Li_in_LiF);
// 	if(!targetLoss_6Li_in_LiF){
// 		ConsoleColorizer::PrintRed("Failed to load TargetEnergyLoss from config file at " + targetEnergyLossPath_6Li_in_LiF + "\n");
// 		return 1;
// 	}

// 	//alphas in LiF target (75 ug/cm^2) for E(alpha) from 0 keV to 20000 keV
// 	std::string targetEnergyLossPath_alpha_in_LiF = "../config/TargetELoss_alpha_in_LiF.conf";
// 	TargetEnergyLoss* targetLoss_alpha_in_LiF = TargetEnergyLoss::LoadFromConfigFile(targetEnergyLossPath_alpha_in_LiF);
// 	if(!targetLoss_alpha_in_LiF){
// 		ConsoleColorizer::PrintRed("Failed to load TargetEnergyLoss from config file at " + targetEnergyLossPath_alpha_in_LiF + "\n");
// 		return 1;
// 	}

// 	//detuerons in LiF target (75 ug/cm^2) for E(deuteron) from 0 keV to 20000 keV
// 	std::string targetenergyLossPath_deuteron_in_LiF = "../config/TargetELoss_deuteron_in_LiF.conf";
// 	TargetEnergyLoss* targetLoss_deuteron_in_LiF = TargetEnergyLoss::LoadFromConfigFile(targetenergyLossPath_deuteron_in_LiF);
// 	if(!targetLoss_deuteron_in_LiF){
// 		ConsoleColorizer::PrintRed("Failed to load TargetEnergyLoss from config file at " + targetenergyLossPath_deuteron_in_LiF + "\n");
// 		return 1;
// 	}

// 	//set up SABRE_DeadLayerModels here:
// 	//6Li in Si (50nm dead layer) for E(6Li) from 0 keV to 10000 keV
// 	std::string deadLayerEnergyLossPath_6Li = "../config/DeadLayerELoss_6Li_in_Si.conf";
// 	SABRE_DeadLayerModel* deadLayerLoss_6Li = SABRE_DeadLayerModel::LoadFromConfigFile(deadLayerEnergyLossPath_6Li);
// 	if(!deadLayerLoss_6Li){
// 		ConsoleColorizer::PrintRed("Failed to load SABRE_DeadLayerModel from config file at " + deadLayerEnergyLossPath_6Li);
// 		return 1;
// 	}

// 	//alphas in Si (50nm dead layer) for E(alpha) from 0 keV to 20000 keV
// 	std::string deadLayerEnergyLossPath_alpha = "../config/DeadLayerELoss_alpha_in_Si.conf";
// 	SABRE_DeadLayerModel* deadLayerLoss_alpha = SABRE_DeadLayerModel::LoadFromConfigFile(deadLayerEnergyLossPath_alpha);
// 	if(!deadLayerLoss_alpha){
// 		ConsoleColorizer::PrintRed("Failed to load SABRE_DeadLayerModel from config file at " + deadLayerEnergyLossPath_alpha + "\n");
// 		return 1;
// 	}

// 	//deuterons in Si (50nm dead layer) for E(deuteron) from 0 keV to 20000 keV
// 	std::string deadLayerEnergyLossPath_deuteron = "../config/DeadLayerELoss_deuteron_in_Si.conf";
// 	SABRE_DeadLayerModel* deadLayerLoss_deuteron = SABRE_DeadLayerModel::LoadFromConfigFile(deadLayerEnergyLossPath_deuteron);
// 	if(!deadLayerLoss_deuteron){
// 		ConsoleColorizer::PrintRed("Failed to load SABRE_DeadLayerModel from config file at " + deadLayerEnergyLossPath_deuteron + "\n");
// 		return 1;
// 	}


// 	//ROOT application for GUI handling:
// 	//TApplication app("SABREsim",&argc,argv);
// 	TH1D* hEnergyLosses = new TH1D("hEnergyLosses","Energy Loss Distribution for 6Li in Si;Energy Loss (keV);Counts",50,0,50);

// 	if(kinX == 2){//kin2mc
// 		std::vector<double> losses;
// 		while(infile >> e1 >> theta1 >> phi1 >> thetacm >> e2 >> theta2 >> phi2){
// 			nevents += 1;
// 			detected1 = false;
// 			detected2 = false;

// 			//beamspot.Spread();//not necessary right now since Set() not working right

// 			outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << thetacm << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << endl;
// 			if(nevents%50000==0) cout << "Processed " << nevents << " events..." << endl;

// 			for(size_t i=0; i<SABRE_Array.size(); i++){
// 				/*///////////////////////////////////////////////
// 				//  check if particle 1 (ejectile) is detected //
// 				///////////////////////////////////////////////*/

// 				//prep variables to hold "smeared" ring and wedge energy after applying resolution
// 				double smearedERing, smearedEWedge;

// 				//get <ring,wedge> pair based on theta,phi
// 				pair<int,int> hit1_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta1*DEG2RAD, phi1*DEG2RAD);
// 				//pair<int,int> hit1_rw = SABRE_Array[i]->GetOffsetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD,{0.000001,0.000001,0.000001});
// 				if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){

// 					//apply target energy loss to e1:
// 					double e1_aftertarget = targetLoss_6Li_in_LiF->ApplyEnergyLoss(e1, theta1);

// 					//apply dead layer energy loss to e1_aftertarget:
// 					Vec3 trajectory;
// 					trajectory.SetVectorSpherical(1,theta1*DEG2RAD,phi1*DEG2RAD);
// 					Vec3 normal = SABRE_Array[i]->GetNormTilted();
// 					normal = normal*(1/normal.Mag());
// 					double e1_afterDeadLayer = deadLayerLoss_6Li->ApplyEnergyLoss(e1_aftertarget, trajectory, normal);

// 					// if(nevents%10000 == 0){
// 					// 	std::cout << "hit1 target_energy_loss = " << abs(e1-e1_aftertarget) << " MeV" << std::endl;
// 					// }

// 					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit1_rw.first,e1_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit1_rw.second,e1_afterDeadLayer,smearedEWedge)){
// 						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
// 						outfile << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
// 						detected1 = true;
// 						hit1+=1;
// 						SABRE_Array_hits[i] += 1;
// 					}
// 				}


// 				/*/////////////////////////////////////////////
// 				//  check if particle 2 (recoil) is detected //
// 				/////////////////////////////////////////////*/

// 				smearedERing = 0.;
// 				smearedEWedge = 0.;

// 				pair<int,int> hit2_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta2*DEG2RAD,phi2*DEG2RAD);
// 				//pair<int,int> hit2_rw = SABRE_Array[i]->GetOffsetTrajectoryRingWedge(theta2*DEG2RAD,phi2*DEG2RAD,{0.000001,0.000001,0.000001});
// 				if(hit2_rw.first != -1 && hit2_rw.second != -1 && !detected2){

// 					//apply target energy loss to e1:
// 					double e2_aftertarget = targetLoss_6Li_in_LiF->ApplyEnergyLoss(e2, theta2);

// 					//apply dead layer energy loss to e1_aftertarget:
// 					Vec3 trajectory;
// 					trajectory.SetVectorSpherical(1,theta2*DEG2RAD,phi2*DEG2RAD);
// 					Vec3 normal = SABRE_Array[i]->GetNormTilted();
// 					normal = normal*(1/normal.Mag());
// 					double e2_afterDeadLayer = deadLayerLoss_6Li->ApplyEnergyLoss(e2_aftertarget, trajectory, normal);

// 					//if(nevents%10000 == 0) std::cout << "hit2 target_energy_loss = " << e2-e2_aftertarget << " MeV" << std::endl;
// 					//hEnergyLosses->Fill(abs(e2-e2_aftertarget)*1000.);//for target energy loss
// 					hEnergyLosses->Fill((e2_afterDeadLayer-e2_aftertarget)*1000.);//for dead layer energy loss

// 					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit2_rw.first,e2_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit2_rw.second,e2_afterDeadLayer,smearedEWedge)){
// 						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit2_rw.first,hit2_rw.second);
// 						outfile << 200+i << "\t" << hit2_rw.first << "\t" << hit2_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
// 						detected2 = true;
// 						hit2+=1;
// 						SABRE_Array_hits[i] += 1;
// 					}
// 				}
// 			}

// 			if(detected1&&detected2){
// 				hitboth += 1;
// 			}
// 			if(detected1&&!detected2){
// 				hit1only += 1;
// 			}
// 			if(!detected1&&detected2){
// 				hit2only += 1;
// 			}

// 			outfile << eoev << endl;

// 		}

// 		TCanvas *c1 = new TCanvas("c1","SABREsim ELoss Dist",800,600);
// 		gStyle->SetOptStat(1111);
// 		hEnergyLosses->Draw();
// 		// app.SetReturnFromRun(true);
// 		// app.Run();
// 		c1->Update();
// 		c1->SaveAs("ELossDist.png");

// 	} else if(kinX == 3){//kin3mc
// 		while(infile >> e1 >> theta1 >> phi1 >> e2 >> theta2 >> phi2 >> e3 >> theta3 >> phi3 >> e4 >> theta4 >> phi4){
// 			nevents += 1;
// 			detected1 = false;
// 			detected3 = false;
// 			detected4 = false;

// 			//outfile << std::format("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",e1,theta1,phi1,e2,theta2,phi2,e3,theta3,phi3,e4,theta4,phi4) << endl;
// 			outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\t" << e3 << "\t" << theta3 << "\t" << phi3 << "\t" << e4 << "\t" << theta4 << "\t" << phi4 << endl;
// 			if(nevents%50000==0) cout << "Processed " << nevents << " events..." << endl;

// 			for(size_t i=0; i<SABRE_Array.size(); i++){
// 				/*///////////////////////////////////////////////
// 				//  check if particle 1 (ejectile) is detected //
// 				///////////////////////////////////////////////*/
// 				double smearedERing, smearedEWedge;

// 				pair<int,int> hit1_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD);
// 				if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){

// 					//apply target energy loss to e1:
// 					double e1_aftertarget = targetLoss_6Li_in_LiF->ApplyEnergyLoss(e1, theta1);

// 					//apply dead layer energy loss to e1_aftertarget:
// 					Vec3 trajectory;
// 					trajectory.SetVectorSpherical(1,theta1*DEG2RAD,phi1*DEG2RAD);
// 					Vec3 normal = SABRE_Array[i]->GetNormTilted();
// 					normal = normal*(1/normal.Mag());
// 					double e1_afterDeadLayer = deadLayerLoss_6Li->ApplyEnergyLoss(e1_aftertarget, trajectory, normal);

// 					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit1_rw.first,e1_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit1_rw.second,e1_afterDeadLayer,smearedEWedge)){
// 						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
// 						//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",100+i,hit1_rw.first,hit1_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
// 						outfile << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
// 						detected1=true;
// 						hit1+=1;
// 						SABRE_Array_hits[i] += 1;
// 					}
// 				}


// 				/*//////////////////////////////////////////
// 				//  check if particle 3 (bu1) is detected //
// 				//////////////////////////////////////////*/
// 				smearedERing = 0.;
// 				smearedEWedge = 0.;

// 				pair<int,int> hit3_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta3*DEG2RAD,phi3*DEG2RAD);

// 				if(hit3_rw.first != -1 && hit3_rw.second != -1 && !detected3){

// 					//apply target energy loss to e3:
// 					double e3_aftertarget = targetLoss_6Li_in_LiF->ApplyEnergyLoss(e3, theta3);

// 					//apply dead layer energy loss to e3_aftertarget:
// 					Vec3 trajectory;
// 					trajectory.SetVectorSpherical(1,theta3*DEG2RAD,phi3*DEG2RAD);
// 					Vec3 normal = SABRE_Array[i]->GetNormTilted();
// 					normal = normal*(1/normal.Mag());
// 					double e3_afterDeadLayer = deadLayerLoss_6Li->ApplyEnergyLoss(e3_aftertarget, trajectory, normal);

// 					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit3_rw.first,e3_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit3_rw.second,e3_afterDeadLayer,smearedEWedge)){
// 						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit3_rw.first,hit3_rw.second);
// 						//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",300+i,hit3_rw.first,hit3_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
// 						outfile << 300+i << "\t" << hit3_rw.first << "\t" << hit3_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
// 						detected3 = true;
// 						hit3 += 1;
// 						SABRE_Array_hits[i] += 1;
// 					}
// 				}


// 				/*//////////////////////////////////////////
// 				//  check if particle 4 (bu2) is detected //
// 				//////////////////////////////////////////*/
// 				smearedERing = 0.;
// 				smearedEWedge = 0.;

// 				pair<int,int> hit4_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta4*DEG2RAD,phi4*DEG2RAD);

// 				if(hit4_rw.first != -1 && hit4_rw.second != -1 && !detected4){

// 					//apply target energy loss to e4:
// 					double e4_aftertarget = targetLoss_6Li_in_LiF->ApplyEnergyLoss(e4, theta4);

// 					//apply dead layer energy loss to e4_aftertarget:
// 					Vec3 trajectory;
// 					trajectory.SetVectorSpherical(1,theta4*DEG2RAD,phi4*DEG2RAD);
// 					Vec3 normal = SABRE_Array[i]->GetNormTilted();
// 					normal = normal*(1/normal.Mag());
// 					double e4_afterDeadLayer = deadLayerLoss_6Li->ApplyEnergyLoss(e4_aftertarget, trajectory, normal);

// 					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit4_rw.first,e4_afterDeadLayer,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit4_rw.second,e4_afterDeadLayer,smearedEWedge)){
// 						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit4_rw.first,hit4_rw.second);
// 						//outfile << std::format("%d\t%d\t%d\t%f\t%f\t%f\t%f",400+i,hit4_rw.first,hit4_rw.second,smearedERing,smearedEWedge,localCoords.GetX(),localCoords.GetY()) << endl;
// 						outfile << 400+i << "\t" << hit4_rw.first << "\t" << hit4_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
// 						detected4 = true;
// 						hit4 += 1;
// 						SABRE_Array_hits[i] += 1;
// 					}
// 				}

// 			}

// 			outfile << eoev << endl;

// 			if(detected3 && detected4){
// 				hit34 += 1;
// 			}
// 			if(detected3&&!detected4){
// 				hitonly3 += 1;
// 			}
// 			if(detected4&&!detected3){
// 				hitonly4 += 1;
// 			}

// 			if((detected1&&!detected3&&!detected4) || (!detected1&&detected3&&!detected4) || (!detected1&&!detected3&&detected4)){
// 				onePartHits += 1;
// 			} else if((detected1&&detected3&&!detected4) || (!detected1&&detected3&&detected4) || (detected1&&!detected3&&detected4)){
// 				twoPartHits += 1;
// 			} else if((detected1&&detected3&&detected4)){
// 				threePartHits += 1;
// 			}
// 		}
// 	} else if(kinX == 4){//kin4mc
// 		while(infile >> e1 >> theta1 >> phi1 >> e2 >> theta2 >> phi2 >> e3 >> theta3 >> phi3 >> e4 >> theta4 >> phi4){
// 			nevents += 1;
// 			detected1 = false;
// 			detected2 = false;
// 			detected3 = false;
// 			detected4 = false;

// 			outfile << e1 << "\t" << theta1 << "\t" << phi1 << "\t" << e2 << "\t" << theta2 << "\t" << phi2 << "\t" << e3 << "\t" << theta3 << "\t" << phi3 << "\t" << e4 << "\t" << theta4 << "\t" << phi4 << endl;
// 			if(nevents%50000==0) cout << "Processed " << nevents << " events..." << endl;

// 			for(size_t i=0; i<SABRE_Array.size(); i++){
// 				double smearedERing, smearedEWedge;

// 				pair<int,int> hit1_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta1*DEG2RAD,phi1*DEG2RAD);
// 				if(hit1_rw.first != -1 && hit1_rw.second != -1 && !detected1){
// 					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit1_rw.first,e1,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit1_rw.second,e1,smearedEWedge)){
// 						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit1_rw.first,hit1_rw.second);
// 						outfile << 100+i << "\t" << hit1_rw.first << "\t" << hit1_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
// 						detected1 = true;
// 						hit1+=1;
// 						SABRE_Array_hits[i] += 1;
// 					}
// 				}

// 				smearedERing = 0.;
// 				smearedEWedge = 0.;
// 				pair<int,int> hit2_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta2*DEG2RAD,phi2*DEG2RAD);
// 				if(hit2_rw.first != -1 && hit2_rw.second != -1 && !detected2){
// 					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit2_rw.first,e2,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit2_rw.second,e2,smearedEWedge)){
// 						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit2_rw.first,hit2_rw.second);	
// 						outfile << 200+i << "\t" << hit2_rw.first << "\t" << hit2_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
// 						detected2 = true;
// 						hit2+=1;
// 						SABRE_Array_hits[i] += 1;
// 					}
// 				}

// 				smearedERing = 0.;
// 				smearedEWedge = 0.;
// 				pair<int,int> hit3_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta3*DEG2RAD,phi3*DEG2RAD);
// 				if(hit3_rw.first != -1 && hit3_rw.second != -1 && !detected3){
// 					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit3_rw.first,e3,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit3_rw.second,e3,smearedEWedge)){
// 						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit3_rw.first,hit3_rw.second);
// 						outfile << 300+i << "\t" << hit3_rw.first << "\t" << hit3_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
// 						detected3 = true;
// 						hit3+=1;
// 						SABRE_Array_hits[i] += 1;
// 					}
// 				}

// 				smearedERing = 0.;
// 				smearedEWedge = 0.;
// 				pair<int,int> hit4_rw = SABRE_Array[i]->GetTrajectoryRingWedge(theta4*DEG2RAD,phi4*DEG2RAD);
// 				if(hit4_rw.first != -1 && hit4_rw.second != -1 && !detected4){
// 					if(SABREARRAY_EnergyResolutionModels[i]->detectEnergyInRing(hit4_rw.first,e4,smearedERing) && SABREARRAY_EnergyResolutionModels[i]->detectEnergyInWedge(hit4_rw.second,e4,smearedEWedge)){
// 						Vec3 localCoords = SABRE_Array[i]->GetHitCoordinatesRandomWiggle(hit4_rw.first,hit4_rw.second);
// 						outfile << 400+i << "\t" << hit4_rw.first << "\t" << hit4_rw.second << "\t" << smearedERing << "\t" << smearedEWedge << "\t" << localCoords.GetX() << "\t" << localCoords.GetY() << endl;
// 						detected4 = true;
// 						hit4+=1;
// 						SABRE_Array_hits[i] += 1;
// 					}
// 				}
// 			}

// 			outfile << eoev << endl;

// 			if((detected1&&!detected2&&!detected3&&!detected4) || (!detected1&&detected2&&!detected3&&!detected4) || (!detected1&&!detected2&&detected3&&!detected4) || (!detected1&&!detected2&&!detected3&&detected4)){
// 				onePartHits += 1;
// 			} else if((detected1&&detected2&&!detected3&&!detected4) || (!detected1&&detected2&&detected3&&!detected4) || (!detected1&&!detected2&&detected3&&detected4) || (detected1&&!detected2&&!detected3&&detected4)){
// 				twoPartHits += 1;
// 			} else if((detected1&&detected2&&detected3&&!detected4) || (!detected1&&detected2&&detected3&&detected4) || (detected1&&!detected2&&detected3&&detected4) || (detected1&&detected2&&!detected3&&detected4)){
// 				threePartHits += 1;
// 			} else if((detected1&&detected2&&detected3&&detected4)){
// 				fourPartHits += 1;
// 			}
// 		}
// 	}

// 	infile.close();
// 	outfile.close();

// 	cout << endl;
// 	cout << "Files are closed. Processed " << nevents << " events.\n" << endl;
// 	if(kinX==2||kinX==3||kinX==4) cout << "Events with ejectile in SABRE: " << hit1 << " (" << float(hit1)*100./float(nevents) << "\%)" << endl;
// 	if(kinX==2) cout << "Events with recoil in SABRE: " << hit2 << " (" << float(hit2)*100./float(nevents) << "\%)\nEvents with only ejectile in SABRE: " << hit1only << " (" << float(hit1only)*100./float(nevents) << "\%)\nEvents with only recoil in SABRE: " << hit2only << " (" << float(hit2only)*100./float(nevents) << "\%)\nEvents with both recoil and ejectile in SABRE: " << hitboth << " (" << float(hitboth)*100./float(nevents) << "\%)" << endl;
// 	if(kinX==3) cout << "Events with (at least) bu1 in SABRE: " << hit3 << "\nEvents with (at least) bu2 in SABRE: " << hit4 << "\nEvents with only bu1 in SABRE: " << hitonly3 << " (" << hitonly3*100./nevents << "\% of all events)" << "\nEvents with only bu2 in SABRE: " << hitonly4 << " (" << hitonly4*100./nevents << "\% of all events)" << "\nEvents with both bu1 and bu2 in SABRE: " << hit34 << " (" << hit34*100./nevents << "\% of all events)" << endl;
// 	if(kinX==4) cout << "Events with bu1 in SABRE: " << hit2 << "\nEvents with bu2 in SABRE: " << hit3 << "\nEvents with bu3 in SABRE: " << hit4 << endl; 
// 	//cout << "Events with both ejectile and recoil in SABRE: " << hitboth << endl << endl;
// 	if(onePartHits>0) cout << "1 Particle Events: " << onePartHits << endl;
// 	if(twoPartHits>0) cout << "2 Particle Events: " << twoPartHits << endl;
// 	if(threePartHits>0) cout << "3 Particle Events: " << threePartHits << endl;
// 	if(fourPartHits>0) cout << "4 Particle Events: " << fourPartHits << endl;
// 	for(int i=0; i<5; i++){
// 		cout << "Detector_" << i << " had total hits = " << SABRE_Array_hits[i] << endl;
// 	}
// 	cout << endl;
// 	cout << "Output can be found here: " << argv[3] << endl;
// 	cout << endl;

// 	for(SABRE_Detector* detector : SABRE_Array){
// 		delete detector;
// 	}
// 	SABRE_Array_hits.clear();

// 	return 1;
// }

int main(int argc, char * argv[]){

	if(argc == 2 && (std::string(argv[1]) == "help" || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")){
		ConsoleColorizer::PrintWhite("Usage:\n\t");
		ConsoleColorizer::PrintYellow("./SABREsim X kinXmcFile.out detOutputFile.det\n\n");
		ConsoleColorizer::PrintWhite("Where ");
		ConsoleColorizer::PrintYellow("X = 2,3,4 ");
		ConsoleColorizer::PrintWhite("for kin2mc, kin3mc, kin4mc input  files\nAnd ");
		ConsoleColorizer::PrintYellow("kinInputFile.out ");
		ConsoleColorizer::PrintWhite("is the filename (or path to) the kinXmc output file\nAnd ");
		ConsoleColorizer::PrintYellow("detOutputFile.det ");
		ConsoleColorizer::PrintWhite("is the filename (or path to) the detection output file (what this code writes to)\n\nExample command:\n");
		ConsoleColorizer::PrintYellow("\t./SABREsim 2 ../kinmc/TEST.out ../det/TEST.det\n\n");
		return 1;
	}

	if(argc!=4){
		ConsoleColorizer::PrintRed("\nError! Please provide input as command line arguments!\nExpected arguments in order of: ");
		ConsoleColorizer::PrintYellow("X kinXmcFile.out detOutputFile.det\n");
		ConsoleColorizer::PrintRed("Where ");
		ConsoleColorizer::PrintYellow("X = 2,3,4 ");
		ConsoleColorizer::PrintRed("for kin2mc, kin3mc, kin4mc input  files\nAnd ");
		ConsoleColorizer::PrintYellow("kinInputFile.out ");
		ConsoleColorizer::PrintRed("is the filename (or path to) the kinXmc output file\nAnd ");
		ConsoleColorizer::PrintYellow("detOutputFile.det ");
		ConsoleColorizer::PrintRed("is the filename (or path to) the detection output file (what this code writes to)\n\nExample command:\n");
		ConsoleColorizer::PrintYellow("\t./SABREsim 2 ../kinmc/TEST.out ../det/TEST.det\n\n");
		return 1;
	}

	int kinX = std::stoi(argv[1]);
	if(kinX != 2 && kinX !=3 && kinX !=4){
		ConsoleColorizer::PrintRed("Error: Invalid kinematics type. Use 2 for 2-body, 3 for 3-body, 4 for 4-body!");
		return 1;
	}

	std::cout << std::endl;
	TString msg = Form("Kin%dmc selected!",kinX);
	ConsoleColorizer::PrintGreen(msg.Data());

	msg = Form("Processing physics data file %s\n", argv[2]);
	ConsoleColorizer::PrintGreen(msg.Data());

	msg = Form("Writing to output file %s\n", argv[3]);
	ConsoleColorizer::PrintGreen(msg.Data());

	try{

		SABREsim sim(kinX, argv[2], argv[3]);
		sim.Run();
		return 0;

	} catch(const std::exception& e){
		ConsoleColorizer::PrintRed(Form("Exception: %s\n",e.what()));
		return 1;
	}
	catch(...){
		ConsoleColorizer::PrintRed("Unknown exception occured during simulation!\n");
		return 1;
	}
}
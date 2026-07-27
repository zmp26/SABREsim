/*

	rays.cpp
	Zach Purcell

	The purpose of this code is to sample across solid angle and toss rays at the SABRE array in order to
	determine the appropriate uncertainties to use on a per-pixel basis.

	A user-defined number of rays are generated evenly across solid angle. These ray trajectories are then
	checked against all SABRE pixels defined and, if intersecting a sensitive portion of SABRE, the "true"
	(read: sampled) theta and phi will be entered into that pixel's 2D histogram of theta vs phi.
	Once all rays have been tossed and checked, the histograms may be fitted with gaussians to yield a best
	fit centroid position and resolution (width sigma). These values can then be used to verify the anglemaps
	and add appropriate angular uncertainties on a per-pixel basis (important for future kinematic fitting).

	Note: This code is compiled and not run as a ROOT macro. Compilation command below:

g++ rays.cpp \
../../src/SABRE_Detector.cpp \
../../src/SABRE_Array.cpp \
-I../../include \
`root-config --cflags --libs` \
-o rays

*/

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <fstream>

#include <TFile.h>
#include <TH2F.h>
#include <TString.h>

#include "../../include/SABRE_Detector.h"
//#include "../../src/SABRE_Detector.cpp"
#include "../../include/SABRE_Array.h"
//#include "../../src/SABRE_Array.cpp"

const std::pair<int,int> offsets[] = {
		{112,40},	//detector0 {ringOffset,wedgeOffset}
		{96,32},	//detector1 {ringOffset,wedgeOffset}
		{80,16},	//detector2 {ringOffset,wedgeOffset}
		{64,24},	//detector3 {ringOffset,wedgeOffset}
		{48,0}		//detector4 {ringOffset,wedgeOffset}
	};

void GenerateAnglesAndUncertainties(){

	long long numThrows = 10000000; //10 million throws!

	double DEG2RAD = M_PI/180.;
	double RAD2DEG = 180./M_PI;

	int nrings = 16;
	int nwedges = 8;

	double INNER_R = 0.0326;
	double OUTER_R = 0.1351;
	double TILT = 40.0;
	double ZDIST = -0.1245;
	double PHI_COVERAGE = 54.4;
	std::vector<double> PHI = {306., 18., 234., 162., 90.};

	SABRE_Array* array = new SABRE_Array();

	for(size_t i=0; i<PHI.size(); i++){
		SABRE_Detector* det = new SABRE_Detector(INNER_R, OUTER_R, PHI_COVERAGE*DEG2RAD, PHI[i]*DEG2RAD, TILT*DEG2RAD, ZDIST);
		array->push_back(det);
	}

	TFile* outfile = new TFile("rays.root", "RECREATE");

	// std::vector<TH2F*> h2_trueAngles((int)PHI.size() * nRings * nWedges, nullptr);
	std::vector<std::vector<TH2F*>> h2_trueAngles(PHI.size());
	std::vector<std::vector<TH2F*>> h2_trueDifferences(PHI.size());
	std::vector<std::vector<std::pair<double,double>>> h2_trueAngles_thetaminmax(PHI.size());
	std::vector<std::vector<std::pair<double,double>>> h2_trueAngles_phiminmax(PHI.size()); 

	for(size_t det=0; det<PHI.size(); det++){
		for(int ring=0; ring<nrings; ring++){
			for(int wedge=0; wedge<nwedges; wedge++){
				TString hname = Form("h2_TrueAngles_Det%zu_Ring%d_Wedge%d", det, ring, wedge);
				TString htitle = Form("True Incidence Angles for Det %zu Ring %d, Wedge %d;#theta_{true} [deg];#phi_{true} [deg]", det, ring, wedge);
				h2_trueAngles[det].push_back(new TH2F(hname, htitle, 360., 90., 180., 1440, 0., 360.));
				hname = Form("h2_TrueDifferences_Det%zu_Ring%d_Wedge%d", det, ring, wedge);
				htitle = Form("True Difference for Det %zu Ring %d Wedge %d;#theta_{pixel} - #theta_{true} [deg];#phi_{pixel} - #phi_{true} [deg]", det, ring, wedge);
				h2_trueDifferences[det].push_back(new TH2F(hname, htitle, 100, -5., 5., 100, -5., 5.));

				h2_trueAngles_thetaminmax[det].push_back({1e9, -1e9});//preload with artificially large and small min max for auto updating
				h2_trueAngles_phiminmax[det].push_back({1e9, -1e9});
			}
		}
	}

	std::random_device rd;
	std::mt19937 gen(rd());

	//we know SABRE covers only ~105deg to ~165deg lab theta, so might as well restrict to this domain for better stats
	double costheta_min = std::cos(100.*DEG2RAD);
	double costheta_max = std::cos(170.*DEG2RAD);
	std::uniform_real_distribution<double> dist_cos_theta(costheta_min, costheta_max);
	std::uniform_real_distribution<double> dist_phi(0., 2.*M_PI);

	std::cout << "Starting Monte Carlo throw loop with uniform solid angle sampling..." << std::endl;
	std::cout << "Will throw " << numThrows << " throws!" << std::endl;

	for(long long i=0; i<numThrows; i++){
		if (i % 1000000 == 0 && i > 0) std::cout << i << " throws processed. (" << i << "/" << numThrows << " = " << i*100./float(numThrows) << "%)" << std::endl;

		double cos_theta = dist_cos_theta(gen);
		double theta_true = std::acos(cos_theta);
		double phi_true = dist_phi(gen);

		int ring = -666;
		int wedge = -666;
		for(size_t det=0; det<PHI.size(); det++){
			std::pair<int, int> hit = array->at(det)->GetTrajectoryRingWedge(theta_true, phi_true);
			if(hit.first >= 0 && hit.first < nrings && hit.second >= 0 && hit.second < nwedges){
				ring = hit.first;
				wedge = hit.second;

				int index = nwedges*ring + wedge;
				h2_trueAngles[det].at(index)->Fill(theta_true*RAD2DEG, phi_true*RAD2DEG);

				auto pixelrw = array->GetDetectorThetaPhi(ring+offsets[det].first, wedge+offsets[det].second);
				if(pixelrw.has_value()){
					double theta_pixel = pixelrw->first;
					double phi_pixel = pixelrw->second;
					h2_trueDifferences[det].at(index)->Fill((theta_pixel - theta_true*RAD2DEG), (phi_pixel - phi_true*RAD2DEG));
					//std::cout << "theta pixel is " << theta_pixel << "\ttheta true is " << theta_true*RAD2DEG << "\n";
					//std::cout << "phi pixel is " << phi_pixel << "\tphi true is " << phi_true*RAD2DEG << "\n\n";
				}

				if(theta_true*RAD2DEG < h2_trueAngles_thetaminmax[det].at(index).first) h2_trueAngles_thetaminmax[det].at(index) = {theta_true*RAD2DEG, h2_trueAngles_thetaminmax[det].at(index).second};
				if(theta_true*RAD2DEG > h2_trueAngles_thetaminmax[det].at(index).second) h2_trueAngles_thetaminmax[det].at(index) = {h2_trueAngles_thetaminmax[det].at(index).first, theta_true*RAD2DEG};

				if(phi_true*RAD2DEG < h2_trueAngles_phiminmax[det].at(index).first) h2_trueAngles_phiminmax[det].at(index) = {phi_true*RAD2DEG, h2_trueAngles_phiminmax[det].at(index).second};
				if(phi_true*RAD2DEG > h2_trueAngles_phiminmax[det].at(index).second) h2_trueAngles_phiminmax[det].at(index) = {h2_trueAngles_phiminmax[det].at(index).first, phi_true*RAD2DEG};
			}

		}
	}

	std::ofstream summaryFile("pixel_uncertainties.txt");
	for(size_t det=0; det<PHI.size(); det++){
		std::string detFilename = "pixel_uncertainties_det" + std::to_string(det) + ".txt";
		std::ofstream detFile(detFilename);
		detFile << "# Global_Ring Global_Wedge Theta_Mean Theta_Sigma Phi_Mean Phi_Sigma\n";

		for(int ring=0; ring<nrings; ring++){
			for(int wedge=0; wedge<nwedges; wedge++){
				int index = nwedges*ring + wedge;

				auto t_mm = h2_trueAngles_thetaminmax[det].at(index);
				auto p_mm = h2_trueAngles_phiminmax[det].at(index);

				// Calculate mean and uniform standard deviation (width / sqrt(12))
				double theta_mean = (t_mm.first + t_mm.second) / 2.0;
				double theta_sigma = (t_mm.second - t_mm.first) / std::sqrt(12.0);

				double phi_mean = (p_mm.first + p_mm.second) / 2.0;
				double phi_sigma = (p_mm.second - p_mm.first) / std::sqrt(12.0);

				// Write to original combined file
				summaryFile << det << " " << ring << " " << wedge << " " 
							<< theta_mean << " " << theta_sigma << " " 
							<< phi_mean << " " << phi_sigma << "\n";

				// Calculate global ring and wedge using offsets for per-detector file
				int global_ring = ring + offsets[det].first;
				int global_wedge = wedge + offsets[det].second;

				// Write to per-detector file (without writing detector index)
				detFile << global_ring << " " << global_wedge << " " 
						<< theta_mean << " " << theta_sigma << " " 
						<< phi_mean << " " << phi_sigma << "\n";
			}
		}
		detFile.close();
	}
	summaryFile.close();

	outfile->Write();
	outfile->Close();

	std::cout << "Finished throwing " << numThrows << " rays at SABRE!" << std::endl;
	std::cout << "Derived pixel means and uncertainties saved to pixel_uncertainties.txt" << std::endl;
}

#ifndef __CINT__
int main() {
	GenerateAnglesAndUncertainties();
	return 0;
}
#endif
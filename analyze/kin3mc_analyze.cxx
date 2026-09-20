#include <fstream>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"

void kin3mc_analyze(const char* filename){

	TMassTable table;
	table.Init("/home/zachpurcell/masstable/masstable.dat");
	
	double massgs_12C = table.GetNuclearMassMeV("C",12);
	double massgs_a = table.GetNuclearMassMeV("He",4);
	double massgs_16O = table.GetNuclearMassMeV("O",16);

	TH1D *hDecay1ThetaCM = new TH1D("hDecay1ThetaCM","ThetaCM of decay particle", 100, 0, TMath::Pi());

	TH1D *hDaughterThetaCM = new TH1D("hDaughterThetaCM", "ThetaCM of daughter particle", 100, 0, TMath::Pi());

	TH1D *hRelativeAngleCM = new TH1D("hRelativeAngleCM", "Relative Angle Between Decay1,Daughter in CM", 100, 0, TMath::Pi());



	ifstream infile(filename);

	double ejE, ejTheta, ejPhi, recE, recTheta, recPhi, bu1E, bu1Theta, bu1Phi, daughterE, daughterTheta, daughterPhi;

	TLorentzVector alpha, carbon, oxygen;

	int count = 0;

	while(infile >> ejE >> ejTheta >> ejPhi >> recE >> recTheta >> recPhi >> bu1E >> bu1Theta >> bu1Phi >> daughterE >> daughterTheta >> daughterPhi){

		double P1 = std::sqrt(bu1E * (bu1E + 2.0*massgs_a));
		double P1x = P1*std::sin(M_PI*bu1Theta/180.)*std::cos(M_PI*bu1Phi/180.);
		double P1y = P1*std::sin(M_PI*bu1Theta/180.)*std::sin(M_PI*bu1Phi/180.);
		double P1z = P1*std::cos(M_PI*bu1Theta/180.);
		alpha.SetPxPyPzE(P1x, P1y, P1z, bu1E+massgs_a);

		double P2 = std::sqrt(bu1E * (bu1E + 2.0*massgs_a));
		double P2x = P1*std::sin(M_PI*bu1Theta/180.)*std::cos(M_PI*bu1Phi/180.);
		double P2y = P1*std::sin(M_PI*bu1Theta/180.)*std::sin(M_PI*bu1Phi/180.);
		double P2z = P1*std::cos(M_PI*bu1Theta/180.);
		carbon.SetPxPyPzE(P2x, P2y, P2z, daughterE+massgs_12C);

		oxygen = alpha + carbon;

		TVector3 boostrecoil = oxygen.BoostVector();

		//boost into CM frame of decay
		alpha.Boost(-boostrecoil);
		carbon.Boost(-boostrecoil);

		//calculate CM angles
		double decay1thetacm = std::acos(P1z/P1);
		double daughterthetacm = std::acos(P2z/P2);
		double relangle = (alpha.Vect()).Angle(carbon.Vect());

		//fill histograms:
		hDecay1ThetaCM->Fill(decay1thetacm);
		hDaughterThetaCM->Fill(daughterthetacm);
		hRelativeAngleCM->Fill(relangle);


		if(count % 50000 == 0){
			std::cout << "Processed " << count << " entries..." << std::endl;
		}

		count += 1;

	}


}
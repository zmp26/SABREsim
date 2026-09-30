#include <fstream>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"

void kin3mc_analyze2(const char* filename){

	TMassTable table;
	table.Init("/home/zachpurcell/masstable/masstable.dat");

	double massgs_6Li = table.GetNuclearMassMeV("Li",6);
	double massgs_a = table.GetNuclearMassMeV("He",4);
	double massgs_d = table.GetNuclearMassMeV("H",2);

	TH2D *hEjectile_ElabVSThetalab = new TH2D("hEjectile_ElabVSThetalab", "hEjectile_ElabVSThetalab", 180, 0, 180, 2100, -1, 20);

	TH2D *hRecoil_ElabVSThetalab = new TH2D("hRecoil_ElabVSThetalab", "hRecoil_ElabVSThetalab", 180, 0, 180, 2100, -1, 20);

	TH2D *hDecay1_ElabVSThetalab = new TH2D("hDecay1_ElabVSThetalab", "hDecay1_ElabVSThetalab", 180, 0, 180, 2100, -1, 20);
	TH2D *hDecay1_ThetalabVSThetaCM = new TH2D("hDecay1_ThetalabVSThetaCM", "hDecay1_ThetalabVSThetaCM;#theta lab(y);#theta CM (x)", 180, 0, 180, 180, 0, 180);

	TH2D *hDaughter_ElabVSThetalab = new TH2D("hDaughter_ElabVSThetalab","hDaughter_ElabVSThetalab", 180, 0, 180, 2100, -1, 20);
	TH2D *hDaughter_ThetalabVSThetaCM = new TH2D("hDaughter_ThetalabVSThetaCM", "hDaughter_ThetalabVSThetaCM;#theta lab(y);#theta CM (x)", 180, 0, 180, 180, 0, 180);

	TH1D *hDecay1ThetaCM = new TH1D("hDecay1ThetaCM","ThetaCM of decay particle", 100, 0, TMath::Pi());
	TH1D *hDecay1CosThetaCM = new TH1D("hDecay1CosThetaCM", "cos(ThetaCM) of decay particle", 100, -1, 1);

	TH1D *hDaughterThetaCM = new TH1D("hDaughterThetaCM", "ThetaCM of daughter particle", 100, 0, TMath::Pi());
	TH1D *hDaughterCosThetaCM = new TH1D("hDaughterCosThetaCM", "cos(ThetaCM) of daughter particle", 100, -1, 1);

	TH1D *hRelativeAngleCM = new TH1D("hRelativeAngleCM", "Relative Angle Between Decay1,Daughter in CM", 100, 0, TMath::Pi());

	ifstream infile(filename);

	double ejE, ejTheta, ejPhi, recE, recTheta, recPhi, bu1E, bu1Theta, bu1Phi, daughterE, daughterTheta, daughterPhi;

	TLorentzVector alpha, deuteron, lithium;

	int count = 0;

	while(infile >> ejE >> ejTheta >> ejPhi >> recE >> recTheta >> recPhi >> bu1E >> bu1Theta >> bu1Phi >> daughterE >> daughterTheta >> daughterPhi){

		hEjectile_ElabVSThetalab->Fill(ejTheta, ejE);
		hRecoil_ElabVSThetalab->Fill(recTheta, recE);
		hDecay1_ElabVSThetalab->Fill(bu1Theta, bu1E);
		hDaughter_ElabVSThetalab->Fill(daughterTheta, daughterE);

		double P1 = std::sqrt(bu1E * (bu1E + 2.0*massgs_a));
		double P1x = P1*std::sin(M_PI*bu1Theta/180.)*std::cos(M_PI*bu1Phi/180.);
		double P1y = P1*std::sin(M_PI*bu1Theta/180.)*std::sin(M_PI*bu1Phi/180.);
		double P1z = P1*std::cos(M_PI*bu1Theta/180.);
		alpha.SetPxPyPzE(P1x, P1y, P1z, bu1E+massgs_a);

		double P2 = std::sqrt(daughterE * (daughterE + 2.0*massgs_d));
		double P2x = P2*std::sin(M_PI*daughterTheta/180.)*std::cos(M_PI*daughterPhi/180.);
		double P2y = P2*std::sin(M_PI*daughterTheta/180.)*std::sin(M_PI*daughterPhi/180.);
		double P2z = P2*std::cos(M_PI*daughterTheta/180.);
		deuteron.SetPxPyPzE(P2x, P2y, P2z, daughterE+massgs_d);

		// lithium = alpha + deuteron;

		// TVector3 boostrecoil = lithium.BoostVector();

		// //boost into CM frame of decay
		// alpha.Boost(-boostrecoil);
		// deuteron.Boost(-boostrecoil);

		// //calculate CM angles
		// //double decay1thetacm = std::acos(P1z/P1);//incorrectly takes dot product with beam axis as z
		// //double daughterthetacm = std::acos(P2z/P2);//incorrectly takes dot product with beam axis as z
		// double decay1thetacm = alpha.Vect().Angle(lithium.Vect());
		// double daughterthetacm = deuteron.Vect().Angle(lithium.Vect());
		// double relangle = (alpha.Vect()).Angle(deuteron.Vect());

		lithium = alpha + deuteron;

		// Parent/recoil direction in the chosen production frame: lab here.
		const TVector3 recoilAxisLab = lithium.Vect().Unit();

		// Boost needed to reach the reconstructed 6Li* rest frame.
		const TVector3 betaLi = lithium.BoostVector();

		// Preserve the lab vectors if you need them later.
		TLorentzVector alphaCM    = alpha;
		TLorentzVector deuteronCM = deuteron;

		// Transform daughters to the 6Li* decay rest frame.
		alphaCM.Boost(-betaLi);
		deuteronCM.Boost(-betaLi);

		// Helicity angles: daughter direction in parent rest frame,
		// referenced to the parent lab-flight direction.
		const double cosThetaAlphaHel = alphaCM.Vect().Unit().Dot(recoilAxisLab);

		const double cosThetaDeuteronHel = deuteronCM.Vect().Unit().Dot(recoilAxisLab);

		const double thetaAlphaHel = std::acos(std::clamp(cosThetaAlphaHel, -1.0, 1.0));

		const double thetaDeuteronHel = std::acos(std::clamp(cosThetaDeuteronHel, -1.0, 1.0));

		// In an exact two-body decay this is pi in the parent rest frame.
		const double relativeAngleCM = alphaCM.Vect().Angle(deuteronCM.Vect());

		//fill histograms:
		hDecay1ThetaCM->Fill(thetaAlphaHel);
		hDecay1CosThetaCM->Fill(cosThetaAlphaHel);
		hDecay1_ThetalabVSThetaCM->Fill(thetaAlphaHel, bu1Theta);

		hDaughterThetaCM->Fill(thetaDeuteronHel);
		hDaughterCosThetaCM->Fill(cosThetaDeuteronHel);
		hDaughter_ThetalabVSThetaCM->Fill(thetaDeuteronHel, daughterTheta);

		hRelativeAngleCM->Fill(relativeAngleCM);


		if(count % 50000 == 0){
			std::cout << "Processed " << count << " entries..." << std::endl;
		}

		count += 1;

	}

}
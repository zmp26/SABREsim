#include <fstream>

#include "TLorentzVector.h"


//set up for 9B->...->paa (kin4mc decay1=p, decay2=a)
void kin4mc_dalitz(const char* filename){

	TMassTable table;
	table.Init("/home/zachpurcell/masstable/masstable.dat");
	double massgs_8Be = table.GetNuclearMassMeV("Be",8);
	double massgs_a = table.GetNuclearMassMeV("He",4);
	double massgs_p = table.GetNuclearMassMeV("H",1);
	double massgs_9B = table.GetNuclearMassMeV("B",9);
	double massgs_5Li = table.GetNuclearMassMeV("Li",5);

	TH2D *hDalitzInvMass = new TH2D("hDalitzInvMass", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
									800, 55.57e6, 55.69e6,
									800, 21.76e6, 21.84e6);

	TH2D *hDalitzInvMass_pa1 = new TH2D("hDalitzInvMass_pa1", "M^{2}_{p+#alpha_1} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha_1}",
									800, 55.57e6, 55.69e6,
									800, 21.76e6, 21.84e6);

	TH2D *hDalitzInvMass_pa2 = new TH2D("hDalitzInvMass_pa2", "M^{2}_{p+#alpha_2} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha_2}",
									800, 55.57e6, 55.69e6,
									800, 21.76e6, 21.84e6);

	TH1D *hInvMass_aa = new TH1D("hInvMass_aa", "hInvMass_aa", 800, std::sqrt(55.57e6), std::sqrt(55.69e6));
	TH1D *hEx8Be = new TH1D("hEx8Be","hEx8Be",175,-1,6);

	TH1D *hInvMass_pa1 = new TH1D("hInvMass_pa1", "hInvMasS_pa1", 800, std::sqrt(21.76e6), std::sqrt(21.84e6));
	TH1D *hEx5Li_pa1 = new TH1D("hEx5Li_pa1","hEx5Li_pa1",175,-1,6);

	TH1D *hInvMass_pa2 = new TH1D("hInvMass_pa2", "hInvMass_pa2", 800, std::sqrt(21.76e6), std::sqrt(21.84e6));
	TH1D *hEx5Li_pa2 = new TH1D("hEx5Li_pa2","hEx5Li_pa2",175,-1,6);

	TH2D *hInvMass_pa1_vs_pa2 = new TH2D("hInvMass_pa1_vs_pa2", "hInvMass_pa1_vs_pa2;pa2;pa1", 800, std::sqrt(21.76e6), std::sqrt(21.84e6), 800, std::sqrt(21.76e6), std::sqrt(21.84e6));

	TH2D *hEx8Be_VS_Ex9B = new TH2D("hEx8Be_VS_Ex9B", "hEx8Be_VS_Ex9B", 1600, -1, 7, 1600, -1, 7);

	ifstream infile(filename);

	double ejE, ejTheta, ejPhi, bu1E, bu1Theta, bu1Phi, bu2E, bu2Theta, bu2Phi, bu3E, bu3Theta, bu3Phi;

	TLorentzVector proton, alpha1, alpha2;
	int count = 0;

	while(infile >> ejE >> ejTheta >> ejPhi >> bu1E >> bu1Theta >> bu1Phi >> bu2E >> bu2Theta >> bu2Phi >> bu3E >> bu3Theta >> bu3Phi){

		//columns 1,2,3 are ejE,ejTheta,ejPhi
		//columns 4,5,6 are bu1E, bu1Theta, bu1Phi
		//columns 7,8,9 are bu2E, bu2Theta, bu2Phi
		//columns 10,11,12 are bu3E, bu3Theta, bu3Phi

		//assume bu1 is proton, but can update this as needed
		double mass_bu1 = massgs_p;
		double Pp = std::sqrt(bu1E * (bu1E + 2.0 * mass_bu1));
		double Pp_x = Pp*std::sin(M_PI*bu1Theta/180.)*std::cos(M_PI*bu1Phi/180.);
		double Pp_y = Pp*std::sin(M_PI*bu1Theta/180.)*std::sin(M_PI*bu1Phi/180.);
		double Pp_z = Pp*std::cos(M_PI*bu1Theta/180.);
		proton.SetPxPyPzE(Pp_x, Pp_y, Pp_z, bu1E+mass_bu1);

		//assume bu2 is alpha, but can update this as needed
		double mass_bu2 = massgs_a;
		double Pa1 = std::sqrt(bu2E * (bu2E + 2.0 * mass_bu2));
		double Pa1_x = Pa1*std::sin(M_PI*bu2Theta/180.)*std::cos(M_PI*bu2Phi/180.);
		double Pa1_y = Pa1*std::sin(M_PI*bu2Theta/180.)*std::sin(M_PI*bu2Phi/180.);
		double Pa1_z = Pa1*std::cos(M_PI*bu2Theta/180.);
		alpha1.SetPxPyPzE(Pa1_x, Pa1_y, Pa1_z, bu2E+mass_bu2);

		//assume bu3 is alpha, but can update this as needed
		double mass_bu3 = massgs_a;
		double Pa2 = std::sqrt(bu3E * (bu3E + 2.0 * mass_bu3));
		double Pa2_x = Pa2*std::sin(M_PI*bu3Theta/180.)*std::cos(M_PI*bu3Phi/180.);
		double Pa2_y = Pa2*std::sin(M_PI*bu3Theta/180.)*std::sin(M_PI*bu3Phi/180.);
		double Pa2_z = Pa2*std::cos(M_PI*bu3Theta/180.);
		alpha2.SetPxPyPzE(Pa2_x, Pa2_y, Pa2_z, bu3E+mass_bu3);


		hInvMass_aa->Fill((alpha1+alpha2).M());
		hEx8Be->Fill((alpha1+alpha2).M() - massgs_8Be);
		hInvMass_pa1->Fill((proton+alpha1).M());
		hEx5Li_pa1->Fill((proton+alpha1).M() - massgs_5Li);
		hInvMass_pa2->Fill((proton+alpha2).M());
		hEx5Li_pa2->Fill((proton+alpha2).M() - massgs_5Li);
		hInvMass_pa1_vs_pa2->Fill((proton+alpha2).M(), (proton+alpha1).M());

		hEx8Be_VS_Ex9B->Fill((proton+alpha1+alpha2).M()-massgs_9B, (alpha1+alpha2).M()-massgs_8Be);

		hDalitzInvMass->Fill((alpha1+alpha2).M2(), (proton+alpha1).M2());
		hDalitzInvMass_pa1->Fill((alpha1+alpha2).M2(), (proton+alpha1).M2());
		hDalitzInvMass->Fill((alpha1+alpha2).M2(), (proton+alpha2).M2());
		hDalitzInvMass_pa2->Fill((alpha1+alpha2).M2(), (proton+alpha2).M2());

		if(count % 100000 == 0){
			std::cout << "Processed " << count << " entries..." << std::endl;
		}

		count += 1;

	}

	infile.close();

	hInvMass_pa1_vs_pa2->Draw("COLZ");

}
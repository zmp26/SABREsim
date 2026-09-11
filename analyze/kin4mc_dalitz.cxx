#include <fstream>

#include "TLorentzVector.h"
#include "TVector3.h"

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

	TH1D *hCosThetaK_A = new TH1D("hCosThetaK_A","hCosThetaK_A", 200, -1, 1);
	TH1D *hCosThetaK_B = new TH1D("hCosThetaK_B","hCosThetaK_B", 200, -1, 1);
	TH2D *hCosThetaK_AvsB = new TH2D("hCosThetaK_AvsB","hCosThetaK_AvsB",200,-1,1,200,-1,1);

	TH1D *hExEt_A = new TH1D("hExEt_A","hExEt_A",100, 0, 1);
	TH1D *hExEt_B = new TH1D("hExEt_B","hExEt_B",100, 0, 1);
	TH2D *hExEt_AvsB = new TH2D("hExEt_AvsB", "hExEt_AvsB", 100, 0, 1, 100, 0, 1);

	TH2D *hExEt_CosThetaK_A = new TH2D("hExEt_CosThetaK_A","hExEt_CosThetaK_A",200,-1,1,100,0,1);
	TH2D *hExEt_CosThetaK_B = new TH2D("hExEt_CosThetaK_B","hExEt_CosThetaK_B",200,-1,1,100,0,1);


	ifstream infile(filename);

	double ejE, ejTheta, ejPhi, bu1E, bu1Theta, bu1Phi, bu2E, bu2Theta, bu2Phi, bu3E, bu3Theta, bu3Phi;

	TLorentzVector proton, alpha1, alpha2;
	TLorentzVector bu1, bu2, bu3;
	int count = 0;
	double ex_A, eT_A, costhetak_A;
	TVector3 k1_A, k2_A, k3_A, kx_A, ky_A;

	double ex_B, eT_B, costhetak_B;
	TVector3 k1_B, k2_B, k3_B, kx_B, ky_B;

	while(infile >> ejE >> ejTheta >> ejPhi >> bu1E >> bu1Theta >> bu1Phi >> bu2E >> bu2Theta >> bu2Phi >> bu3E >> bu3Theta >> bu3Phi){

		//for kin4mc format apa:
		//double mass_bu1 = massgs_a;
		//double mass_bu2 = massgs_p;
		//double mass_bu3 = massgs_a;

		//for kin4mc format paa:
		double mass_bu1 = massgs_p;
		double mass_bu2 = massgs_a;
		double mass_bu3 = massgs_a;

		//columns 1,2,3 are ejE,ejTheta,ejPhi
		//columns 4,5,6 are bu1E, bu1Theta, bu1Phi
		//columns 7,8,9 are bu2E, bu2Theta, bu2Phi
		//columns 10,11,12 are bu3E, bu3Theta, bu3Phi

		double P1 = std::sqrt(bu1E * (bu1E + 2.0 * mass_bu1));
		double P1x = P1*std::sin(M_PI*bu1Theta/180.)*std::cos(M_PI*bu1Phi/180.);
		double P1y = P1*std::sin(M_PI*bu1Theta/180.)*std::sin(M_PI*bu1Phi/180.);
		double P1z = P1*std::cos(M_PI*bu1Theta/180.);
		bu1.SetPxPyPzE(P1x, P1y, P1z, bu1E+mass_bu1);

		double P2 = std::sqrt(bu2E * (bu2E + 2.0 * mass_bu2));
		double P2x = P2*std::sin(M_PI*bu2Theta/180.)*std::cos(M_PI*bu2Phi/180.);
		double P2y = P2*std::sin(M_PI*bu2Theta/180.)*std::sin(M_PI*bu2Phi/180.);
		double P2z = P2*std::cos(M_PI*bu2Theta/180.);
		bu2.SetPxPyPzE(P2x, P2y, P2z, bu2E+mass_bu2);

		double P3 = std::sqrt(bu3E * (bu3E + 2.0 * mass_bu3));
		double P3x = P3*std::sin(M_PI*bu3Theta/180.)*std::cos(M_PI*bu3Phi/180.);
		double P3y = P3*std::sin(M_PI*bu3Theta/180.)*std::sin(M_PI*bu3Phi/180.);
		double P3z = P3*std::cos(M_PI*bu3Theta/180.);
		bu3.SetPxPyPzE(P3x, P3y, P3z, bu3E+mass_bu3);

		//decay1 = p (p+8Be or dem if decay1=p)
		proton = bu1;
		alpha1 = bu2;
		alpha2 = bu3;

		//decay1 = a (a+5Li or dem if decay1=a)
		// alpha1 = bu1;
		// proton = bu2;
		// alpha2 = bu3;

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


		//prep variables for jacobi calculations:
		TLorentzVector recoil = proton + alpha1 + alpha2;
		TVector3 boostrecoil = recoil.BoostVector();
		proton.Boost(-boostrecoil); //boost into CM frame of 9B recoil
		alpha1.Boost(-boostrecoil); //boost into CM frame of 9B recoil
		alpha2.Boost(-boostrecoil); //boost into CM frame of 9B recoil


		//jacobi coordinates toy model:
		//T system, take first alpha column to be alpha1 and second alpha column to be alpha2
		k1_A = alpha1.Vect();
		k2_A = alpha2.Vect();
		k3_A = proton.Vect();

		// kx_A = (massgs_a*k1_A - massgs_a*k2_A)*(1./(massgs_a+massgs_a));
		// ky_A = (massgs_p*(k1_A+k2_A) - (massgs_a+massgs_a)*k3_A)*(1./(massgs_a+massgs_a+massgs_p));
		kx_A = 0.5*(k1_A - k2_A);
		ky_A = -k3_A;

		ex_A = (massgs_a+massgs_a)*kx_A.Mag2()/(2*massgs_a*massgs_a);
		eT_A = (alpha1.E() - massgs_a + alpha2.E() - massgs_a + proton.E() - massgs_p);

		costhetak_A = (kx_A.Dot(ky_A))/(kx_A.Mag()*ky_A.Mag());

		hCosThetaK_A->Fill(costhetak_A);
		hExEt_A->Fill(ex_A/eT_A);
		hExEt_CosThetaK_A->Fill(costhetak_A, ex_A/eT_A);


		//T system, take second alpha column to be alpha1 and first alpha column to be alpha2
		k1_B = alpha2.Vect();
		k2_B = alpha1.Vect();
		k3_B = proton.Vect();

		//kx_B = (massgs_a*k1_B - massgs_a*k2_B)*(1./(massgs_a+massgs_a));
		//ky_B = (massgs_p*(k1_A+k2_A) - (massgs_a+massgs_a)*k3_A)*(1./(massgs_a+massgs_a+massgs_p));
		kx_B = 0.5*(k1_B - k2_B);
		ky_B = -k3_B;

		ex_B = (massgs_a+massgs_a)*kx_B.Mag2()*(1./(2*massgs_a*massgs_a));
		eT_B = (alpha1.E() - massgs_a + alpha2.E() - massgs_a + proton.E() - massgs_p);

		costhetak_B = (kx_B.Dot(ky_B))/(kx_B.Mag()*ky_B.Mag());

		hCosThetaK_B->Fill(costhetak_B);
		hExEt_B->Fill(ex_B/eT_B);
		hExEt_CosThetaK_B->Fill(costhetak_B, ex_B/eT_B);

		hCosThetaK_AvsB->Fill(costhetak_B, costhetak_A);
		hExEt_AvsB->Fill(ex_B/eT_B, ex_A/eT_A);

		// std::cout << "k1_A = (" << k1_A[0] << ", " << k1_A[1] << ", " << k1_A[2] << ")\n";
		// std::cout << "k1_B = (" << k1_B[0] << ", " << k1_B[1] << ", " << k1_B[2] << ")\n";
		// std::cout << "k1_A - k1_B = (" << k1_A[0] - k1_B[0] << ", " << k1_A[1] - k1_B[1] << ", " << k1_A[2] - k1_B[2] << ")\n";
		// std::cout << "|k1_A - k1_B| = " << (k1_A-k1_B).Mag() << "\n";
		// std::cout << "k2_A = (" << k2_A[0] << ", " << k2_A[1] << ", " << k2_A[2] << ")\n";
		// std::cout << "k2_B = (" << k2_B[0] << ", " << k2_B[1] << ", " << k2_B[2] << ")\n";
		// std::cout << "k2_A - k2_B = (" << k2_A[0] - k2_B[0] << ", " << k2_A[1] - k2_B[1] << ", " << k2_A[2] - k2_B[2] << ")\n";
		// std::cout << "|k2_A - k2_B| = " << (k2_A-k2_B).Mag() << "\n\n";

		if(count % 100000 == 0){
			std::cout << "Processed " << count << " entries..." << std::endl;
		}

		count += 1;

	}

	infile.close();

	hInvMass_pa1_vs_pa2->Draw("COLZ");

}
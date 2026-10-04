#include <fstream>
#include <string>

#include "TLorentzVector.h"
#include "TVector3.h"

//helper function for jacobi kx, ky coordinates given k1, k2, k3, and masses m1, m2, m3:
std::pair<TVector3,TVector3> getJacobiMomenta(TVector3 k1, TVector3 k2, TVector3 k3, double m1, double m2, double m3){

	TVector3 kx, ky;

	kx = (m2*k1 - m1*k2)*(1./(m1+m2));
	ky = (m3*(k1+k2) - (m1+m2)*k3)*(1./(m1+m2+m3));

	return std::make_pair(kx,ky);
}

//set up for 9B->...->paa (kin4mc decay1=p, decay2=a)
//channel
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

		//for kin4mc format paa:
		// double mass_bu1 = massgs_p;
		// double mass_bu2 = massgs_a;
		// double mass_bu3 = massgs_a;

		//for kin4mc format apa:
		double mass_bu1 = massgs_a;
		double mass_bu2 = massgs_p;
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
		// proton = bu1;
		// alpha1 = bu2;
		// alpha2 = bu3;

		//decay1 = a (a+5Li or dem if decay1=a)
		alpha1 = bu1;
		proton = bu2;
		alpha2 = bu3;

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

//set up for 9B->...->paa (kin4mc decay1=p, decay2=a)
//channel = "paa" -> decay1 = p, decay2 = a (9B -> p+8Be, 8Be->a+a)
//channel = "apa" -> decay1 = a, decay2 = p (9B -> a+5Li, 5Li->p+a)
void kin4mc_dalitz(const char* filename, const char* outfilename, const std::string& channel = "paa"){

	TMassTable table;
	table.Init("/home/zachpurcell/masstable/masstable.dat");
	double massgs_8Be = table.GetNuclearMassMeV("Be",8);
	double massgs_a = table.GetNuclearMassMeV("He",4);
	double massgs_p = table.GetNuclearMassMeV("H",1);
	double massgs_9B = table.GetNuclearMassMeV("B",9);
	double massgs_5Li = table.GetNuclearMassMeV("Li",5);

	TFile *outfile = new TFile(outfilename, "RECREATE");


	TH2D *hDalitzInvMass = new TH2D("hDalitzInvMass", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
									800, 55.57e6, 55.69e6,
									800, 21.76e6, 21.84e6);
	hDalitzInvMass->SetDirectory(outfile);

	TH2D *hDalitzInvMass_pa1 = new TH2D("hDalitzInvMass_pa1", "M^{2}_{p+#alpha_1} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha_1}",
									800, 55.57e6, 55.69e6,
									800, 21.76e6, 21.84e6);
	hDalitzInvMass_pa1->SetDirectory(outfile);

	TH2D *hDalitzInvMass_pa2 = new TH2D("hDalitzInvMass_pa2", "M^{2}_{p+#alpha_2} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha_2}",
									800, 55.57e6, 55.69e6,
									800, 21.76e6, 21.84e6);
	hDalitzInvMass_pa2->SetDirectory(outfile);

	TH1D *hInvMass_aa = new TH1D("hInvMass_aa", "hInvMass_aa", 800, std::sqrt(55.57e6), std::sqrt(55.69e6));
	hInvMass_aa->SetDirectory(outfile);
	TH1D *hEx8Be = new TH1D("hEx8Be","hEx8Be",175,-1,6);
	hEx8Be->SetDirectory(outfile);

	TH1D *hInvMass_pa1 = new TH1D("hInvMass_pa1", "hInvMasS_pa1", 800, std::sqrt(21.76e6), std::sqrt(21.84e6));
	hInvMass_pa1->SetDirectory(outfile);
	TH1D *hEx5Li_pa1 = new TH1D("hEx5Li_pa1","hEx5Li_pa1",175,-1,6);
	hEx5Li_pa1->SetDirectory(outfile);

	TH1D *hInvMass_pa2 = new TH1D("hInvMass_pa2", "hInvMass_pa2", 800, std::sqrt(21.76e6), std::sqrt(21.84e6));
	hInvMass_pa2->SetDirectory(outfile);
	TH1D *hEx5Li_pa2 = new TH1D("hEx5Li_pa2","hEx5Li_pa2",175,-1,6);
	hEx5Li_pa2->SetDirectory(outfile);

	TH2D *hInvMass_pa1_vs_pa2 = new TH2D("hInvMass_pa1_vs_pa2", "hInvMass_pa1_vs_pa2;pa2;pa1", 800, std::sqrt(21.76e6), std::sqrt(21.84e6), 800, std::sqrt(21.76e6), std::sqrt(21.84e6));
	hInvMass_pa1_vs_pa2->SetDirectory(outfile);

	TH2D *hEx8Be_VS_Ex9B = new TH2D("hEx8Be_VS_Ex9B", "hEx8Be_VS_Ex9B", 1600, -1, 7, 1600, -1, 7);
	hEx8Be_VS_Ex9B->SetDirectory(outfile);

	//jacobi T histograms

	TH1D *hCosThetaK_T = new TH1D("hCosThetaK_T","hCosThetaK_T", 200, -1, 1);
	hCosThetaK_T->SetDirectory(outfile);

	TH1D *hExEt_T = new TH1D("hExEt_T","hExEt_T",100, 0, 1);
	hExEt_T->SetDirectory(outfile);

	TH1D *hEyEt_T = new TH1D("hEyEt_T", "hEyEt_T", 100, 0, 1);
	hEyEt_T->SetDirectory(outfile);

	TH2D *hExEt_CosThetaK_T = new TH2D("hExEt_CosThetaK_T","hExEt_CosThetaK_T",200,-1,1,100,0,1);
	hExEt_CosThetaK_T->SetDirectory(outfile);

	TH2D *hEyEt_CosThetaK_T = new TH2D("hEyEt_CosThetaK_T","hEyEt_CosThetaK_T",200,-1,1,100,0,1);
	hEyEt_CosThetaK_T->SetDirectory(outfile);

	TH2D *hExEt_Et_T = new TH2D("hExEt_Et_T", "hExEt_Et_T", 1000, 0, 10, 100, 0, 1);
	hExEt_Et_T->SetDirectory(outfile);

	TH2D *hExSPS_ExEt_T = new TH2D("hExSPS_ExEt_T","hExSPS_ExEt_T",100,0,1,1400,0,7);
	hExSPS_ExEt_T->SetDirectory(outfile);

	//jacobi Y1 histograms

	TH1D *hCosThetaK_Y1 = new TH1D("hCosThetaK_Y1","hCosThetaK_Y1", 200, -1, 1);
	hCosThetaK_Y1->SetDirectory(outfile);

	TH1D *hExEt_Y1 = new TH1D("hExEt_Y1","hExEt_Y1",100, 0, 1);
	hExEt_Y1->SetDirectory(outfile);

	TH1D *hEyEt_Y1 = new TH1D("hEyEt_Y1", "hEyEt_Y1", 100, 0, 1);
	hEyEt_Y1->SetDirectory(outfile);

	TH2D *hExEt_CosThetaK_Y1 = new TH2D("hExEt_CosThetaK_Y1","hExEt_CosThetaK_Y1",200,-1,1,100,0,1);
	hExEt_CosThetaK_Y1->SetDirectory(outfile);

	TH2D *hEyEt_CosThetaK_Y1 = new TH2D("hEyEt_CosThetaK_Y1","hEyEt_CosThetaK_Y1",200,-1,1,100,0,1);
	hEyEt_CosThetaK_Y1->SetDirectory(outfile);

	TH2D *hExEt_Et_Y1 = new TH2D("hExEt_Et_Y1", "hExEt_Et_Y1", 1000, 0, 10, 100, 0, 1);
	hExEt_Et_Y1->SetDirectory(outfile);

	TH2D *hExSPS_ExEt_Y1 = new TH2D("hExSPS_ExEt_Y1","hExSPS_ExEt_Y1",100,0,1,1400,0,7);
	hExSPS_ExEt_Y1->SetDirectory(outfile);

	//jacobi Y2 histograms

	TH1D *hCosThetaK_Y2 = new TH1D("hCosThetaK_Y2","hCosThetaK_Y2", 200, -1, 1);
	hCosThetaK_Y2->SetDirectory(outfile);

	TH1D *hExEt_Y2 = new TH1D("hExEt_Y2","hExEt_Y2",100, 0, 1);
	hExEt_Y2->SetDirectory(outfile);

	TH1D *hEyEt_Y2 = new TH1D("hEyEt_Y2", "hEyEt_Y2", 100, 0, 1);
	hEyEt_Y2->SetDirectory(outfile);

	TH2D *hExEt_CosThetaK_Y2 = new TH2D("hExEt_CosThetaK_Y2","hExEt_CosThetaK_Y2",200,-1,1,100,0,1);
	hExEt_CosThetaK_Y2->SetDirectory(outfile);

	TH2D *hEyEt_CosThetaK_Y2 = new TH2D("hEyEt_CosThetaK_Y2","hEyEt_CosThetaK_Y2",200,-1,1,100,0,1);
	hEyEt_CosThetaK_Y2->SetDirectory(outfile);

	TH2D *hExEt_Et_Y2 = new TH2D("hExEt_Et_Y2", "hExEt_Et_Y2", 1000, 0, 10, 100, 0, 1);
	hExEt_Et_Y2->SetDirectory(outfile);

	TH2D *hExSPS_ExEt_Y2 = new TH2D("hExSPS_ExEt_Y2","hExSPS_ExEt_Y2",100,0,1,1400,0,7);
	hExSPS_ExEt_Y2->SetDirectory(outfile);

	//jacobi Y histograms

	TH1D *hCosThetaK_Y = new TH1D("hCosThetaK_Y","hCosThetaK_Y", 200, -1, 1);
	hCosThetaK_Y->SetDirectory(outfile);

	TH1D *hExEt_Y = new TH1D("hExEt_Y","hExEt_Y",100, 0, 1);
	hExEt_Y->SetDirectory(outfile);

	TH1D *hEyEt_Y = new TH1D("hEyEt_Y", "hEyEt_Y", 100, 0, 1);
	hEyEt_Y->SetDirectory(outfile);

	TH2D *hExEt_CosThetaK_Y = new TH2D("hExEt_CosThetaK_Y","hExEt_CosThetaK_Y",200,-1,1,100,0,1);
	hExEt_CosThetaK_Y->SetDirectory(outfile);

	TH2D *hEyEt_CosThetaK_Y = new TH2D("hEyEt_CosThetaK_Y","hEyEt_CosThetaK_Y",200,-1,1,100,0,1);
	hEyEt_CosThetaK_Y->SetDirectory(outfile);

	TH2D *hExEt_Et_Y = new TH2D("hExEt_Et_Y", "hExEt_Et_Y", 1000, 0, 10, 100, 0, 1);
	hExEt_Et_Y->SetDirectory(outfile);

	TH2D *hExSPS_ExEt_Y = new TH2D("hExSPS_ExEt_Y","hExSPS_ExEt_Y",100,0,1,1400,0,7);
	hExSPS_ExEt_Y->SetDirectory(outfile);

	//CM angles (CM of recoil break up, so 9B CM)
	TH1D *hProtonThetaCM = new TH1D("hProtonThetaCM","hProtonThetaCM", 90, 0, 180);
	hProtonThetaCM->SetDirectory(outfile);

	TH1D *hProtonPhiCM = new TH1D("hProtonPhiCM", "hProtonPhiCM", 180, 0, 360);
	hProtonPhiCM->SetDirectory(outfile);

	TH1D *hAlphaThetaCM = new TH1D("hAlphaThetaCM", "hAlphaThetaCM", 90, 0, 180);
	hAlphaThetaCM->SetDirectory(outfile);

	TH1D *hAlphaPhiCM = new TH1D("hAlphaPhiCM", "hAlphaPhiCM", 180, 0, 360);
	hAlphaPhiCM->SetDirectory(outfile);


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

		double mass_bu1, mass_bu2, mass_bu3;

		if(channel == "paa"){
			mass_bu1 = massgs_p;
			mass_bu2 = massgs_a;
			mass_bu3 = massgs_a;
		} else if(channel == "apa"){
			mass_bu1 = massgs_a;
			mass_bu2 = massgs_p;
			mass_bu3 = massgs_a;
		} else {
			std::cout << "Invalid channel!" << std::endl;
			return;
		}

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

		if(channel == "paa"){//decay1=p, decay2=a (9B->p+8Be, 8Be->a+a)
			proton = bu1;
			alpha1 = bu2;
			alpha2 = bu3;
		} else if(channel == "apa"){//decay1=a, decay2=p (9B->a+5Li, 5Li->p+a)
			alpha1 = bu1;
			proton = bu2;
			alpha2 = bu3;
		}

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

		double RecEx = recoil.M() - massgs_9B;


		//jacobi T system with 1=alpha1, 2=alpha2, 3=proton
		double mu_x = massgs_a*massgs_a/(2.*massgs_a);
		double mu_y = massgs_p*2.*massgs_a/(massgs_p+massgs_a+massgs_a);

		TVector3 k1_T, k2_T, k3_T;
		k1_T = alpha1.Vect();
		k2_T = alpha2.Vect();
		k3_T = proton.Vect();

		double protonThetaCM = std::acos(proton.Vect()[2]/proton.Vect().Mag())*180./M_PI;
		double protonPhiCM = std::atan2(proton.Vect()[1],proton.Vect()[0])*180./M_PI;
		if(protonPhiCM < 0) protonPhiCM += 360.;
		hProtonThetaCM->Fill(protonThetaCM);
		hProtonPhiCM->Fill(protonPhiCM);

		double alpha1ThetaCM = std::acos(alpha1.Vect()[2]/alpha1.Vect().Mag())*180./M_PI;
		double alpha1PhiCM = std::atan2(alpha1.Vect()[1],alpha1.Vect()[0])*180./M_PI;
		if(alpha1PhiCM < 0) alpha1PhiCM += 360.;
		hAlphaThetaCM->Fill(alpha1ThetaCM);
		hAlphaPhiCM->Fill(alpha1PhiCM);

		double alpha2ThetaCM = std::acos(alpha2.Vect()[2]/alpha2.Vect().Mag())*180./M_PI;
		double alpha2PhiCM = std::atan2(alpha2.Vect()[1],alpha2.Vect()[0])*180./M_PI;
		if(alpha2PhiCM < 0) alpha2PhiCM += 360.;
		hAlphaThetaCM->Fill(alpha2ThetaCM);
		hAlphaPhiCM->Fill(alpha2PhiCM);


		std::pair<TVector3,TVector3> jacobiMomenta = getJacobiMomenta(k1_T, k2_T, k3_T, massgs_a, massgs_a, massgs_p);
		TVector3 kx_T = jacobiMomenta.first;
		TVector3 ky_T = jacobiMomenta.second;

		double kxTMag = kx_T.Mag();
		double kyTMag = ky_T.Mag();
		double E_xT = kxTMag*kxTMag/(2.*mu_x); //(mass_a+mass_a)*kx_A.Mag2()/(2.*mass_a*mass_a);
		double E_yT = kyTMag*kyTMag/(2.*mu_y); //(mass_p+mass_a+mass_a)*kyMag*kyMag/(2*(mass_p*2*mass_a));
		double E_Txy_T = E_xT + E_yT;	
		double costhetakT = (kx_T.Dot(ky_T))/(kx_T.Mag()*ky_T.Mag());

		hCosThetaK_T->Fill(costhetakT);
		hCosThetaK_T->Fill(-costhetakT);//costhetak is negated under transformation of 1,2 -> 2,1, must fill both
		hExEt_T->Fill(E_xT/E_Txy_T);
		hEyEt_T->Fill(E_yT/E_Txy_T);
		hExEt_Et_T->Fill(E_Txy_T, E_xT/E_Txy_T);
		hExEt_CosThetaK_T->Fill(costhetakT, E_xT/E_Txy_T);
		hExEt_CosThetaK_T->Fill(-costhetakT, E_xT/E_Txy_T);//costhetak is negated under transformation of 1,2 -> 2,1, must fill both
		hEyEt_CosThetaK_T->Fill(costhetakT, E_yT/E_Txy_T);
		hEyEt_CosThetaK_T->Fill(-costhetakT, E_yT/E_Txy_T);//costhetak is negated under transformation of 1,2 -> 2,1, must fill both
		hExSPS_ExEt_T->Fill(E_xT/E_Txy_T, RecEx);

		//jacobi Y system, with 1=alpha1, 2=proton, 3=alpha2
		mu_x = (massgs_a*massgs_p)/(massgs_a + massgs_p);
		mu_y = (massgs_a*(massgs_a+massgs_p))/(massgs_a+massgs_p+massgs_a);

		TVector3 k1_Y, k2_Y, k3_Y;
		k1_Y = alpha1.Vect();
		k2_Y = proton.Vect();
		k3_Y = alpha2.Vect();

		jacobiMomenta = getJacobiMomenta(k1_Y, k2_Y, k3_Y, massgs_a, massgs_p, massgs_a);
		TVector3 kx_Y1 = jacobiMomenta.first;
		TVector3 ky_Y1 = jacobiMomenta.second;

		double kxY1Mag = kx_Y1.Mag();
		double kyY1Mag = ky_Y1.Mag();
		double E_xY1 = kxY1Mag*kxY1Mag/(2.*mu_x);
		double E_yY1 = kyY1Mag*kyY1Mag/(2.*mu_y);
		double E_Txy_Y1 = E_xY1 + E_yY1;
		double costhetakY1 = (kx_Y1.Dot(ky_Y1))/(kxY1Mag*kyY1Mag);

		hCosThetaK_Y1->Fill(costhetakY1);//y1
		hCosThetaK_Y->Fill(costhetakY1);//cumulative
		hExEt_Y1->Fill(E_xY1/E_Txy_Y1);//y1
		hExEt_Y->Fill(E_xY1/E_Txy_Y1);//cumulative
		hEyEt_Y1->Fill(E_yY1/E_Txy_Y1);//y1
		hEyEt_Y->Fill(E_yY1/E_Txy_Y1);//cumulative
		hExEt_Et_Y1->Fill(E_Txy_Y1, E_xY1/E_Txy_Y1);
		hExEt_Et_Y->Fill(E_Txy_Y1, E_xY1/E_Txy_Y1);
		hExEt_CosThetaK_Y1->Fill(costhetakY1, E_xY1/E_Txy_Y1);//y1
		hExEt_CosThetaK_Y->Fill(costhetakY1, E_xY1/E_Txy_Y1);//cumulative
		hEyEt_CosThetaK_Y1->Fill(costhetakY1, E_yY1/E_Txy_Y1);//y1
		hEyEt_CosThetaK_Y->Fill(costhetakY1, E_yY1/E_Txy_Y1);//cumulative
		hExSPS_ExEt_Y1->Fill(E_xY1/E_Txy_Y1, RecEx);
		hExSPS_ExEt_Y->Fill(E_xY1/E_Txy_Y1, RecEx);

		//jacobi Y system, with 1 = alpha2, 2 = proton, 3 = alpha1
		//mux, muy do not change here but calculate anyway
		mu_x = (massgs_a*massgs_p)/(massgs_a+massgs_p);
		mu_y = (massgs_a*(massgs_a+massgs_p))/(massgs_a+massgs_p+massgs_a);

		k1_Y = alpha2.Vect();
		k2_Y = proton.Vect();
		k3_Y = alpha1.Vect();

		jacobiMomenta = getJacobiMomenta(k1_Y, k2_Y, k3_Y, massgs_a, massgs_p, massgs_a);
		TVector3 kx_Y2 = jacobiMomenta.first;
		TVector3 ky_Y2 = jacobiMomenta.second;

		double kxY2Mag = kx_Y2.Mag();
		double kyY2Mag = ky_Y2.Mag();
		double E_xY2 = kxY2Mag*kxY2Mag/(2.*mu_x);
		double E_yY2 = kyY2Mag*kyY2Mag/(2.*mu_y);
		double E_Txy_Y2 = E_xY2 + E_yY2;
		double costhetakY2 = (kx_Y2.Dot(ky_Y2))/(kxY2Mag*kyY2Mag);

		hCosThetaK_Y2->Fill(costhetakY2);//y2
		hCosThetaK_Y->Fill(costhetakY2);//cumulative
		hExEt_Y2->Fill(E_xY2/E_Txy_Y2);//y2
		hExEt_Y->Fill(E_xY2/E_Txy_Y2);//cumulative
		hEyEt_Y2->Fill(E_yY2/E_Txy_Y2);//y2
		hEyEt_Y->Fill(E_yY2/E_Txy_Y2);//cumulative
		hExEt_Et_Y2->Fill(E_Txy_Y2, E_xY2/E_Txy_Y2);
		hExEt_Et_Y->Fill(E_Txy_Y2, E_xY2/E_Txy_Y2);
		hExEt_CosThetaK_Y2->Fill(costhetakY2, E_xY2/E_Txy_Y2);//y2
		hExEt_CosThetaK_Y->Fill(costhetakY2, E_xY2/E_Txy_Y2);//cumulative
		hEyEt_CosThetaK_Y2->Fill(costhetakY2, E_yY2/E_Txy_Y2);//y2
		hEyEt_CosThetaK_Y->Fill(costhetakY2, E_yY2/E_Txy_Y2);//cumulative
		hExSPS_ExEt_Y2->Fill(E_xY2/E_Txy_Y2, RecEx);
		hExSPS_ExEt_Y->Fill(E_xY2/E_Txy_Y2, RecEx);

		//get helicty angle here:
		//	parent frame: 			rest frame of recoil (recoil = p + alpha1 + alpha2)
		//	intermediate: 			alpha1 + alpha2
		//
		//calculation follows "General Properties of Three-body Decays" by Curtis Meyer, Carnegie Mellon University
		//the paper uses the direction of the lorentz boost that takes parent CM frame -> a1+a2 rest frame.
		//Thus, we can say the reference direction is  -beta_aa_in_parent_CM
		TLorentzVector protonCM = proton;
		TLorentzVector alpha1CM = alpha1;
		TLorentzVector alpha2CM = alpha2;

		// TVector3 betaLabToParentCM = -recoil.BoostVector();
		// protonCM.Boost(betaLabToParentCM);
		// alpha1CM.Boost(betaLabToParentCM);
		// alpha2CM.Boost(betaLabToParentCM);

		TLorentzVector intermediateCM = alpha1CM + alpha2CM;
		TVector3 betaCMtoIntermediateCM = -intermediateCM.BoostVector();

		double cosThetaH_a1 = -666.;
		double cosThetaH_a2 = -666.;
		double thetaH_a1_deg = -666.;
		double thetaH_a2_deg = -666.;

		if(betaCMtoIntermediateCM.Mag() > 1e-12 && alpha1CM.P() > 1e-12 && alpha2CM.P() > 1e-12 ){
			//TVector3 helicityAxis = betaCMtoIntermediateCM.Unit();
			TVector3 helicityAxis = intermediateCM.Vect().Unit();
			TLorentzVector alpha1AA = alpha1CM;
			TLorentzVector alpha2AA = alpha2CM;

			alpha1AA.Boost(betaCMtoIntermediateCM);
			alpha2AA.Boost(betaCMtoIntermediateCM);

			cosThetaH_a1 = std::max(-1.,std::min(1., alpha1AA.Vect().Unit().Dot(helicityAxis)));
			cosThetaH_a2 = std::max(-1.,std::min(1., alpha2AA.Vect().Unit().Dot(helicityAxis)));

			thetaH_a1_deg = std::acos(cosThetaH_a1)*180./M_PI;
			thetaH_a2_deg = std::acos(cosThetaH_a2)*180./M_PI;
		}


		if(count % 100000 == 0){
			std::cout << "Processed " << count << " entries..." << std::endl;
		}

		count += 1;

	}

	infile.close();

	//hInvMass_pa1_vs_pa2->Draw("COLZ");
	outfile->Write();
	outfile->Close();
}
//purpose of this code: analyze data from 3par and 2par 10B(3He,a)9B experiment (or sim) that
//						has already been analyzed by SABREanalyze_Mult3.cxx (B10ha_SABREPID and B10ha_SABREPID_Mult2)

//expects the PID_Mult3 or PID_Mult2 tree to exist in the file


#include <utility>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TDirectory.h"

//helper function to check if value is within bounds
bool IsInside(double value, std::pair<double,double> bounds){
	if(bounds.first > bounds.second) return false;
	return (value >= bounds.first && value <= bounds.second);
}

void DalitzRegions(TTree *tree);

//declare channel regions here:
//p+8Be(0+)
std::pair<double, double> p8Be0plus_sig1 = {55.5746e6, 55.5750e6};
std::pair<double, double> p8Be0plus_sig2 = {55.5744e6, 55.5752e6};
std::pair<double, double> p8Be0plus_sig3 = {55.5742e6, 55.5754e6};
bool p8Be0plus_sig1_check = false, p8Be0plus_sig2_check = false, p8Be0plus_sig3_check = false;

//p+8Be(2+)
std::pair<double, double> p8Be2plus_sig1 = {55.610e6, 55.630e6};
std::pair<double, double> p8Be2plus_sig2 = {55.600e6, 55.639e6};
std::pair<double, double> p8Be2plus_sig3 = {55.591e6, 55.649e6};
bool p8Be2plus_sig1_check = false, p8Be2plus_sig2_check = false, p8Be2plus_sig3_check = false;

//a+5Li
std::pair<double, double> a5Ligs_sig1 = {21.7818e6, 21.7915e6};
std::pair<double, double> a5Ligs_sig2 = {21.7769e6, 21.7964e6};
std::pair<double, double> a5Ligs_sig3 = {21.7721e6, 21.8013e6};
bool a5Ligs_sig1_check = false, a5Ligs_sig2_check = false, a5Ligs_sig3_check = false;

TTree *intree = nullptr;
TFile *infile = nullptr;

//use this function to analyze from a file that is not open in your ROOT session
//filename = string containing relative path to file
//SABREmult = 2 or 3 depending on mult2/mult3 data (changes expected tree name)
void DalitzRegions(TString filename, int SABREmult){
	
	if(SABREmult != 2 && SABREmult != 3){
		std::cout << "wrong SABREmult, try 2 or 3" << std::endl;
		return;
	}

	infile = TFile::Open(filename, "READ");
	if(!infile || infile->IsZombie()){
		std::cout << "bad file" << std::endl;
		return;
	}

	intree = (TTree*)infile->Get(SABREmult == 3 ? "PID_Mult3" : "PID_Mult2");
	if(!intree){
		std::cout << "bad tree in file" << std::endl;
		infile->Close();
		return;
	}
	//intree->SetDirectory(0);

	//infile->Close();

	DalitzRegions(intree);
}

void DalitzRegions(TTree *tree){
	if(!tree){
		std::cout << "passed tree is bad" << std::endl;
		return;
	}

	long nentries = tree->GetEntries();

	std::cout << "tree has " << nentries << " entries\n";

	//prepare output file:
	TString outfilename = Form("%s_regions.root", tree->GetName());
	TFile *outfile = new TFile(outfilename,"RECREATE");
	if(!outfile || outfile->IsZombie()){
		std::cout << "failed to create output file" << std::endl;
		return;
	}

	//tree->SetDirectory(outfile);

	//prepare histograms:
	//all channels//////////////////////////////////////////////////////////////////////////////////////////////////
	TH2D *hDalitzInvMass = new TH2D("hDalitzInvMass", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
									800, 55.57e6, 55.69e6,
									800, 21.76e6, 21.84e6);

	//SABREsum histos
	TH2D *hSABREsumE_vs_ExSPS = new TH2D("hSABREsumE_vs_ExSPS","hSABREsumE_vs_ExSPS", 1400, 0, 7, 1400, 0, 7);

	//9B Ex histos:
	TH1D *hExSPS = new TH1D("hExSPS", "hExSPS", 1600, -1, 7);
	TH1D *hRecEx = new TH1D("hRecEx", "hRecEx", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS = new TH2D("hRecEx_VS_ExSPS", "hRecEx_VS_ExSPS", 1600, -1, 7, 1600, -1, 7);

	//jacobi histos
		TH1D *hCosThetaK_T = new TH1D("hCosThetaK_T","hCosThetaK_T", 200, -1, 1);
		TH1D *hExEt_T = new TH1D("hExEt_T","hExEt_T",100, 0, 1);
		TH1D *hEyEt_T = new TH1D("hEyEt_T", "hEyEt_T", 100, 0, 1);
		TH2D *hExEt_CosThetaK_T = new TH2D("hExEt_CosThetaK_T","hExEt_CosThetaK_T",200,-1,1,100,0,1);
		TH2D *hEyEt_CosThetaK_T = new TH2D("hEyEt_CosThetaK_T","hEyEt_CosThetaK_T",200,-1,1,100,0,1);
		TH2D *hExEt_Et_T = new TH2D("hExEt_Et_T", "hExEt_Et_T", 1000, 0, 10, 100, 0, 1);
		TH2D *hExSPS_ExEt_T = new TH2D("hExSPS_ExEt_T","hExSPS_ExEt_T",100,0,1,1400,0,7);

		TH1D *hCosThetaK_Y1 = new TH1D("hCosThetaK_Y1","hCosThetaK_Y1", 200, -1, 1);
		TH1D *hExEt_Y1 = new TH1D("hExEt_Y1","hExEt_Y1",100, 0, 1);
		TH1D *hEyEt_Y1 = new TH1D("hEyEt_Y1", "hEyEt_Y1", 100, 0, 1);
		TH2D *hExEt_CosThetaK_Y1 = new TH2D("hExEt_CosThetaK_Y1","hExEt_CosThetaK_Y1",200,-1,1,100,0,1);
		TH2D *hEyEt_CosThetaK_Y1 = new TH2D("hEyEt_CosThetaK_Y1","hEyEt_CosThetaK_Y1",200,-1,1,100,0,1);
		TH2D *hExEt_Et_Y1 = new TH2D("hExEt_Et_Y1", "hExEt_Et_Y1", 1000, 0, 10, 100, 0, 1);
		TH2D *hExSPS_ExEt_Y1 = new TH2D("hExSPS_ExEt_Y1","hExSPS_ExEt_Y1",100,0,1,1400,0,7);

		TH1D *hCosThetaK_Y2 = new TH1D("hCosThetaK_Y2","hCosThetaK_Y2", 200, -1, 1);
		TH1D *hExEt_Y2 = new TH1D("hExEt_Y2","hExEt_Y2",100, 0, 1);
		TH1D *hEyEt_Y2 = new TH1D("hEyEt_Y2", "hEyEt_Y2", 100, 0, 1);
		TH2D *hExEt_CosThetaK_Y2 = new TH2D("hExEt_CosThetaK_Y2","hExEt_CosThetaK_Y2",200,-1,1,100,0,1);
		TH2D *hEyEt_CosThetaK_Y2 = new TH2D("hEyEt_CosThetaK_Y2","hEyEt_CosThetaK_Y2",200,-1,1,100,0,1);
		TH2D *hExEt_Et_Y2 = new TH2D("hExEt_Et_Y2", "hExEt_Et_Y2", 1000, 0, 10, 100, 0, 1);
		TH2D *hExSPS_ExEt_Y2 = new TH2D("hExSPS_ExEt_Y2","hExSPS_ExEt_Y2",100,0,1,1400,0,7);

	TH2D *h2_a5Ligs_VS_p8Be0plus = new TH2D("h2_a5Ligs_VS_p8Be0plus", "h2_a5Ligs_VS_p8Be0plus", 3, 0.5, 3.5, 3, 0.5, 3.5);
	h2_a5Ligs_VS_p8Be0plus->GetXaxis()->SetTitle("p+8Be(0^{+}) sigmas");
	h2_a5Ligs_VS_p8Be0plus->GetYaxis()->SetTitle("a+5Li(gs) sigmas");
	h2_a5Ligs_VS_p8Be0plus->SetDirectory(outfile);

	// TH2D *h2_a5Ligs_VS_p8Be0plus_STRICT = new TH2D("h2_a5Ligs_VS_p8Be0plus_STRICT", "h2_a5Ligs_VS_p8Be0plus_STRICT", 3, 0.5, 3.5, 3, 0.5, 3.5);
	// h2_a5Ligs_VS_p8Be0plus_STRICT->GetXaxis()->SetTitle("p+8Be(0^{+}) sigmas");
	// h2_a5Ligs_VS_p8Be0plus_STRICT->GetYaxis()->SetTitle("a+5Li(gs) sigmas");
	// h2_a5Ligs_VS_p8Be0plus_STRICT->SetDirectory(outfile);
	
	TH2D *h2_a5Ligs_VS_p8Be2plus = new TH2D("h2_a5Ligs_VS_p8Be2plus", "h2_a5Ligs_VS_p8Be2plus", 3, 0.5, 3.5, 3, 0.5, 3.5);
	h2_a5Ligs_VS_p8Be2plus->GetXaxis()->SetTitle("p+8Be(2^{+}) sigmas");
	h2_a5Ligs_VS_p8Be2plus->GetYaxis()->SetTitle("a+5Li(gs) sigmas");
	h2_a5Ligs_VS_p8Be2plus->SetDirectory(outfile);

	// TH2D *h2_a5Ligs_VS_p8Be2plus_STRICT = new TH2D("h2_a5Ligs_VS_p8Be2plus_STRICT", "h2_a5Ligs_VS_p8Be2plus_STRICT", 3, 0.5, 3.5, 3, 0.5, 3.5);
	// h2_a5Ligs_VS_p8Be2plus_STRICT->GetXaxis()->SetTitle("p+8Be(2^{+}) sigmas");
	// h2_a5Ligs_VS_p8Be2plus_STRICT->GetYaxis()->SetTitle("a+5Li(gs) sigmas");
	// h2_a5Ligs_VS_p8Be2plus_STRICT->SetDirectory(outfile);

	// TH2D *h2_p8Be0plus_VS_p8Be2plus = new TH2D("h2_p8Be0plus_VS_p8Be2plus", "h2_p8Be0plus_VS_p8Be2plus", 3, 0.5, 3.5, 3, 0.5, 3.5);
	// h2_p8Be0plus_VS_p8Be2plus->GetXaxis()->SetTitle("p+8Be(2^{+}) sigmas");
	// h2_p8Be0plus_VS_p8Be2plus->GetYaxis()->SetTitle("a+5Li(gs) sigmas");
	// h2_p8Be0plus_VS_p8Be2plus->SetDirectory(outfile);


	hDalitzInvMass->SetDirectory(outfile);

	hSABREsumE_vs_ExSPS->SetDirectory(outfile);

	hExSPS->SetDirectory(outfile);
	hRecEx->SetDirectory(outfile);
	hRecEx_VS_ExSPS->SetDirectory(outfile);

	hCosThetaK_T->SetDirectory(outfile);
	hExEt_T->SetDirectory(outfile);
	hEyEt_T->SetDirectory(outfile);
	hExEt_CosThetaK_T->SetDirectory(outfile);
	hEyEt_CosThetaK_T->SetDirectory(outfile);
	hExEt_Et_T->SetDirectory(outfile);
	hExSPS_ExEt_T->SetDirectory(outfile);

	hCosThetaK_Y1->SetDirectory(outfile);
	hExEt_Y1->SetDirectory(outfile);
	hEyEt_Y1->SetDirectory(outfile);
	hExEt_CosThetaK_Y1->SetDirectory(outfile);
	hEyEt_CosThetaK_Y1->SetDirectory(outfile);
	hExEt_Et_Y1->SetDirectory(outfile);
	hExSPS_ExEt_Y1->SetDirectory(outfile);

	hCosThetaK_Y2->SetDirectory(outfile);
	hExEt_Y2->SetDirectory(outfile);
	hEyEt_Y2->SetDirectory(outfile);
	hExEt_CosThetaK_Y2->SetDirectory(outfile);
	hEyEt_CosThetaK_Y2->SetDirectory(outfile);
	hExEt_Et_Y2->SetDirectory(outfile);
	hExSPS_ExEt_Y2->SetDirectory(outfile);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//p+8Be(0+)////////////////////////////////////////////////////////////////////////////////////////////////////

	//dalitz plots
	TH2D *hDalitzInvMass_p8Be0_1 = new TH2D("hDalitzInvMass_p8Be0_1", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);
	TH2D *hDalitzInvMass_p8Be0_2 = new TH2D("hDalitzInvMass_p8Be0_2", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);
	TH2D *hDalitzInvMass_p8Be0_3 = new TH2D("hDalitzInvMass_p8Be0_3", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);

	//SABREsum histos
	TH2D *hSABREsumE_vs_ExSPS_p8Be0_1 = new TH2D("hSABREsumE_vs_ExSPS_p8Be0_1","hSABREsumE_vs_ExSPS_p8Be0_1", 1400, 0, 7, 1400, 0, 7);
	TH2D *hSABREsumE_vs_ExSPS_p8Be0_2 = new TH2D("hSABREsumE_vs_ExSPS_p8Be0_2","hSABREsumE_vs_ExSPS_p8Be0_2", 1400, 0, 7, 1400, 0, 7);
	TH2D *hSABREsumE_vs_ExSPS_p8Be0_3 = new TH2D("hSABREsumE_vs_ExSPS_p8Be0_3","hSABREsumE_vs_ExSPS_p8Be0_3", 1400, 0, 7, 1400, 0, 7);

	//9B Ex histos:
	TH1D *hExSPS_p8Be0_1 = new TH1D("hExSPS_p8Be0_1", "hExSPS_p8Be0_1", 1600, -1, 7);
	TH1D *hRecEx_p8Be0_1 = new TH1D("hRecEx_p8Be0_1", "hRecEx_p8Be0_1", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_p8Be0_1 = new TH2D("hRecEx_VS_ExSPS_p8Be0_1", "hRecEx_VS_ExSPS_p8Be0_1", 1600, -1, 7, 1600, -1, 7);

	TH1D *hExSPS_p8Be0_2 = new TH1D("hExSPS_p8Be0_2", "hExSPS_p8Be0_2", 1600, -1, 7);
	TH1D *hRecEx_p8Be0_2 = new TH1D("hRecEx_p8Be0_2", "hRecEx_p8Be0_2", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_p8Be0_2 = new TH2D("hRecEx_VS_ExSPS_p8Be0_2", "hRecEx_VS_ExSPS_p8Be0_2", 1600, -1, 7, 1600, -1, 7);

	TH1D *hExSPS_p8Be0_3 = new TH1D("hExSPS_p8Be0_3", "hExSPS_p8Be0_3", 1600, -1, 7);
	TH1D *hRecEx_p8Be0_3 = new TH1D("hRecEx_p8Be0_3", "hRecEx_p8Be0_3", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_p8Be0_3 = new TH2D("hRecEx_VS_ExSPS_p8Be0_3", "hRecEx_VS_ExSPS_p8Be0_3", 1600, -1, 7, 1600, -1, 7);

	//jacobi histos
		//1sig
			//T
				TH1D *hCosThetaK_T_p8Be0_1 = new TH1D("hCosThetaK_T_p8Be0_1","hCosThetaK_T_p8Be0_1", 200, -1, 1);
				TH1D *hExEt_T_p8Be0_1 = new TH1D("hExEt_T_p8Be0_1","hExEt_T_p8Be0_1",100, 0, 1);
				TH1D *hEyEt_T_p8Be0_1 = new TH1D("hEyEt_T_p8Be0_1", "hEyEt_T_p8Be0_1", 100, 0, 1);
				TH2D *hExEt_CosThetaK_T_p8Be0_1 = new TH2D("hExEt_CosThetaK_T_p8Be0_1","hExEt_CosThetaK_T_p8Be0_1",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_T_p8Be0_1 = new TH2D("hEyEt_CosThetaK_T_p8Be0_1","hEyEt_CosThetaK_T_p8Be0_1",200,-1,1,100,0,1);
				TH2D *hExEt_Et_T_p8Be0_1 = new TH2D("hExEt_Et_T_p8Be0_1", "hExEt_Et_T_p8Be0_1", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_T_p8Be0_1 = new TH2D("hExSPS_ExEt_T_p8Be0_1","hExSPS_ExEt_T_p8Be0_1",100,0,1,1400,0,7);
			//Y1
				TH1D *hCosThetaK_Y1_p8Be0_1 = new TH1D("hCosThetaK_Y1_p8Be0_1","hCosThetaK_Y1_p8Be0_1", 200, -1, 1);
				TH1D *hExEt_Y1_p8Be0_1 = new TH1D("hExEt_Y1_p8Be0_1","hExEt_Y1_p8Be0_1",100, 0, 1);
				TH1D *hEyEt_Y1_p8Be0_1 = new TH1D("hEyEt_Y1_p8Be0_1", "hEyEt_Y1_p8Be0_1", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y1_p8Be0_1 = new TH2D("hExEt_CosThetaK_Y1_p8Be0_1","hExEt_CosThetaK_Y1_p8Be0_1",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y1_p8Be0_1 = new TH2D("hEyEt_CosThetaK_Y1_p8Be0_1","hEyEt_CosThetaK_Y1_p8Be0_1",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y1_p8Be0_1 = new TH2D("hExEt_Et_Y1_p8Be0_1", "hExEt_Et_Y1_p8Be0_1", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y1_p8Be0_1 = new TH2D("hExSPS_ExEt_Y1_p8Be0_1","hExSPS_ExEt_Y1_p8Be0_1",100,0,1,1400,0,7);
			//Y2
				TH1D *hCosThetaK_Y2_p8Be0_1 = new TH1D("hCosThetaK_Y2_p8Be0_1","hCosThetaK_Y2_p8Be0_1", 200, -1, 1);
				TH1D *hExEt_Y2_p8Be0_1 = new TH1D("hExEt_Y2_p8Be0_1","hExEt_Y2_p8Be0_1",100, 0, 1);
				TH1D *hEyEt_Y2_p8Be0_1 = new TH1D("hEyEt_Y2_p8Be0_1", "hEyEt_Y2_p8Be0_1", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y2_p8Be0_1 = new TH2D("hExEt_CosThetaK_Y2_p8Be0_1","hExEt_CosThetaK_Y2_p8Be0_1",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y2_p8Be0_1 = new TH2D("hEyEt_CosThetaK_Y2_p8Be0_1","hEyEt_CosThetaK_Y2_p8Be0_1",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y2_p8Be0_1 = new TH2D("hExEt_Et_Y2_p8Be0_1", "hExEt_Et_Y2_p8Be0_1", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y2_p8Be0_1 = new TH2D("hExSPS_ExEt_Y2_p8Be0_1","hExSPS_ExEt_Y2_p8Be0_1",100,0,1,1400,0,7);
		//2sig
			//T
				TH1D *hCosThetaK_T_p8Be0_2 = new TH1D("hCosThetaK_T_p8Be0_2","hCosThetaK_T_p8Be0_2", 200, -1, 1);
				TH1D *hExEt_T_p8Be0_2 = new TH1D("hExEt_T_p8Be0_2","hExEt_T_p8Be0_2",100, 0, 1);
				TH1D *hEyEt_T_p8Be0_2 = new TH1D("hEyEt_T_p8Be0_2", "hEyEt_T_p8Be0_2", 100, 0, 1);
				TH2D *hExEt_CosThetaK_T_p8Be0_2 = new TH2D("hExEt_CosThetaK_T_p8Be0_2","hExEt_CosThetaK_T_p8Be0_2",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_T_p8Be0_2 = new TH2D("hEyEt_CosThetaK_T_p8Be0_2","hEyEt_CosThetaK_T_p8Be0_2",200,-1,1,100,0,1);
				TH2D *hExEt_Et_T_p8Be0_2 = new TH2D("hExEt_Et_T_p8Be0_2", "hExEt_Et_T_p8Be0_2", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_T_p8Be0_2 = new TH2D("hExSPS_ExEt_T_p8Be0_2","hExSPS_ExEt_T_p8Be0_2",100,0,1,1400,0,7);
			//Y1
				TH1D *hCosThetaK_Y1_p8Be0_2 = new TH1D("hCosThetaK_Y1_p8Be0_2","hCosThetaK_Y1_p8Be0_2", 200, -1, 1);
				TH1D *hExEt_Y1_p8Be0_2 = new TH1D("hExEt_Y1_p8Be0_2","hExEt_Y1_p8Be0_2",100, 0, 1);
				TH1D *hEyEt_Y1_p8Be0_2 = new TH1D("hEyEt_Y1_p8Be0_2", "hEyEt_Y1_p8Be0_2", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y1_p8Be0_2 = new TH2D("hExEt_CosThetaK_Y1_p8Be0_2","hExEt_CosThetaK_Y1_p8Be0_2",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y1_p8Be0_2 = new TH2D("hEyEt_CosThetaK_Y1_p8Be0_2","hEyEt_CosThetaK_Y1_p8Be0_2",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y1_p8Be0_2 = new TH2D("hExEt_Et_Y1_p8Be0_2", "hExEt_Et_Y1_p8Be0_2", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y1_p8Be0_2 = new TH2D("hExSPS_ExEt_Y1_p8Be0_2","hExSPS_ExEt_Y1_p8Be0_2",100,0,1,1400,0,7);
			//Y2
				TH1D *hCosThetaK_Y2_p8Be0_2 = new TH1D("hCosThetaK_Y2_p8Be0_2","hCosThetaK_Y2_p8Be0_2", 200, -1, 1);
				TH1D *hExEt_Y2_p8Be0_2 = new TH1D("hExEt_Y2_p8Be0_2","hExEt_Y2_p8Be0_2",100, 0, 1);
				TH1D *hEyEt_Y2_p8Be0_2 = new TH1D("hEyEt_Y2_p8Be0_2", "hEyEt_Y2_p8Be0_2", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y2_p8Be0_2 = new TH2D("hExEt_CosThetaK_Y2_p8Be0_2","hExEt_CosThetaK_Y2_p8Be0_2",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y2_p8Be0_2 = new TH2D("hEyEt_CosThetaK_Y2_p8Be0_2","hEyEt_CosThetaK_Y2_p8Be0_2",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y2_p8Be0_2 = new TH2D("hExEt_Et_Y2_p8Be0_2", "hExEt_Et_Y2_p8Be0_2", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y2_p8Be0_2 = new TH2D("hExSPS_ExEt_Y2_p8Be0_2","hExSPS_ExEt_Y2_p8Be0_2",100,0,1,1400,0,7);
		//3sig
			//T
				TH1D *hCosThetaK_T_p8Be0_3 = new TH1D("hCosThetaK_T_p8Be0_3","hCosThetaK_T_p8Be0_3", 200, -1, 1);
				TH1D *hExEt_T_p8Be0_3 = new TH1D("hExEt_T_p8Be0_3","hExEt_T_p8Be0_3",100, 0, 1);
				TH1D *hEyEt_T_p8Be0_3 = new TH1D("hEyEt_T_p8Be0_3", "hEyEt_T_p8Be0_3", 100, 0, 1);
				TH2D *hExEt_CosThetaK_T_p8Be0_3 = new TH2D("hExEt_CosThetaK_T_p8Be0_3","hExEt_CosThetaK_T_p8Be0_3",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_T_p8Be0_3 = new TH2D("hEyEt_CosThetaK_T_p8Be0_3","hEyEt_CosThetaK_T_p8Be0_3",200,-1,1,100,0,1);
				TH2D *hExEt_Et_T_p8Be0_3 = new TH2D("hExEt_Et_T_p8Be0_3", "hExEt_Et_T_p8Be0_3", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_T_p8Be0_3 = new TH2D("hExSPS_ExEt_T_p8Be0_3","hExSPS_ExEt_T_p8Be0_3",100,0,1,1400,0,7);
			//Y1
				TH1D *hCosThetaK_Y1_p8Be0_3 = new TH1D("hCosThetaK_Y1_p8Be0_3","hCosThetaK_Y1_p8Be0_3", 200, -1, 1);
				TH1D *hExEt_Y1_p8Be0_3 = new TH1D("hExEt_Y1_p8Be0_3","hExEt_Y1_p8Be0_3",100, 0, 1);
				TH1D *hEyEt_Y1_p8Be0_3 = new TH1D("hEyEt_Y1_p8Be0_3", "hEyEt_Y1_p8Be0_3", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y1_p8Be0_3 = new TH2D("hExEt_CosThetaK_Y1_p8Be0_3","hExEt_CosThetaK_Y1_p8Be0_3",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y1_p8Be0_3 = new TH2D("hEyEt_CosThetaK_Y1_p8Be0_3","hEyEt_CosThetaK_Y1_p8Be0_3",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y1_p8Be0_3 = new TH2D("hExEt_Et_Y1_p8Be0_3", "hExEt_Et_Y1_p8Be0_3", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y1_p8Be0_3 = new TH2D("hExSPS_ExEt_Y1_p8Be0_3","hExSPS_ExEt_Y1_p8Be0_3",100,0,1,1400,0,7);
			//Y2
				TH1D *hCosThetaK_Y2_p8Be0_3 = new TH1D("hCosThetaK_Y2_p8Be0_3","hCosThetaK_Y2_p8Be0_3", 200, -1, 1);
				TH1D *hExEt_Y2_p8Be0_3 = new TH1D("hExEt_Y2_p8Be0_3","hExEt_Y2_p8Be0_3",100, 0, 1);
				TH1D *hEyEt_Y2_p8Be0_3 = new TH1D("hEyEt_Y2_p8Be0_3", "hEyEt_Y2_p8Be0_3", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y2_p8Be0_3 = new TH2D("hExEt_CosThetaK_Y2_p8Be0_3","hExEt_CosThetaK_Y2_p8Be0_3",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y2_p8Be0_3 = new TH2D("hEyEt_CosThetaK_Y2_p8Be0_3","hEyEt_CosThetaK_Y2_p8Be0_3",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y2_p8Be0_3 = new TH2D("hExEt_Et_Y2_p8Be0_3", "hExEt_Et_Y2_p8Be0_3", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y2_p8Be0_3 = new TH2D("hExSPS_ExEt_Y2_p8Be0_3","hExSPS_ExEt_Y2_p8Be0_3",100,0,1,1400,0,7);

	//prep directories:
	TDirectory *p8Be0plus = outfile->mkdir("p8Be0plus");
	TDirectory *p8Be0plus_1 = p8Be0plus->mkdir("p8Be0_1sig");
	TDirectory *p8Be0plus_2 = p8Be0plus->mkdir("p8Be0_2sig");
	TDirectory *p8Be0plus_3 = p8Be0plus->mkdir("p8Be0_3sig");

	//set histograms owner to proper directory
	hDalitzInvMass_p8Be0_1->SetDirectory(p8Be0plus_1);
	hDalitzInvMass_p8Be0_2->SetDirectory(p8Be0plus_2);
	hDalitzInvMass_p8Be0_3->SetDirectory(p8Be0plus_3);

	hSABREsumE_vs_ExSPS_p8Be0_1->SetDirectory(p8Be0plus_1);
	hSABREsumE_vs_ExSPS_p8Be0_2->SetDirectory(p8Be0plus_2);
	hSABREsumE_vs_ExSPS_p8Be0_3->SetDirectory(p8Be0plus_3);

	hExSPS_p8Be0_1->SetDirectory(p8Be0plus_1);
	hRecEx_p8Be0_1->SetDirectory(p8Be0plus_1);
	hRecEx_VS_ExSPS_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExSPS_p8Be0_2->SetDirectory(p8Be0plus_2);
	hRecEx_p8Be0_2->SetDirectory(p8Be0plus_2);
	hRecEx_VS_ExSPS_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExSPS_p8Be0_3->SetDirectory(p8Be0plus_3);
	hRecEx_p8Be0_3->SetDirectory(p8Be0plus_3);
	hRecEx_VS_ExSPS_p8Be0_3->SetDirectory(p8Be0plus_3);

	hCosThetaK_T_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExEt_T_p8Be0_1->SetDirectory(p8Be0plus_1);
	hEyEt_T_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExEt_CosThetaK_T_p8Be0_1->SetDirectory(p8Be0plus_1);
	hEyEt_CosThetaK_T_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExEt_Et_T_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExSPS_ExEt_T_p8Be0_1->SetDirectory(p8Be0plus_1);

	hCosThetaK_T_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExEt_T_p8Be0_2->SetDirectory(p8Be0plus_2);
	hEyEt_T_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExEt_CosThetaK_T_p8Be0_2->SetDirectory(p8Be0plus_2);
	hEyEt_CosThetaK_T_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExEt_Et_T_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExSPS_ExEt_T_p8Be0_2->SetDirectory(p8Be0plus_2);

	hCosThetaK_T_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExEt_T_p8Be0_3->SetDirectory(p8Be0plus_3);
	hEyEt_T_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExEt_CosThetaK_T_p8Be0_3->SetDirectory(p8Be0plus_3);
	hEyEt_CosThetaK_T_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExEt_Et_T_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExSPS_ExEt_T_p8Be0_3->SetDirectory(p8Be0plus_3);

	hCosThetaK_Y1_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExEt_Y1_p8Be0_1->SetDirectory(p8Be0plus_1);
	hEyEt_Y1_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExEt_CosThetaK_Y1_p8Be0_1->SetDirectory(p8Be0plus_1);
	hEyEt_CosThetaK_Y1_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExEt_Et_Y1_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExSPS_ExEt_Y1_p8Be0_1->SetDirectory(p8Be0plus_1);

	hCosThetaK_Y1_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExEt_Y1_p8Be0_2->SetDirectory(p8Be0plus_2);
	hEyEt_Y1_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExEt_CosThetaK_Y1_p8Be0_2->SetDirectory(p8Be0plus_2);
	hEyEt_CosThetaK_Y1_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExEt_Et_Y1_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExSPS_ExEt_Y1_p8Be0_2->SetDirectory(p8Be0plus_2);

	hCosThetaK_Y1_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExEt_Y1_p8Be0_3->SetDirectory(p8Be0plus_3);
	hEyEt_Y1_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExEt_CosThetaK_Y1_p8Be0_3->SetDirectory(p8Be0plus_3);
	hEyEt_CosThetaK_Y1_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExEt_Et_Y1_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExSPS_ExEt_Y1_p8Be0_3->SetDirectory(p8Be0plus_3);

	hCosThetaK_Y2_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExEt_Y2_p8Be0_1->SetDirectory(p8Be0plus_1);
	hEyEt_Y2_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExEt_CosThetaK_Y2_p8Be0_1->SetDirectory(p8Be0plus_1);
	hEyEt_CosThetaK_Y2_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExEt_Et_Y2_p8Be0_1->SetDirectory(p8Be0plus_1);
	hExSPS_ExEt_Y2_p8Be0_1->SetDirectory(p8Be0plus_1);

	hCosThetaK_Y2_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExEt_Y2_p8Be0_2->SetDirectory(p8Be0plus_2);
	hEyEt_Y2_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExEt_CosThetaK_Y2_p8Be0_2->SetDirectory(p8Be0plus_2);
	hEyEt_CosThetaK_Y2_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExEt_Et_Y2_p8Be0_2->SetDirectory(p8Be0plus_2);
	hExSPS_ExEt_Y2_p8Be0_2->SetDirectory(p8Be0plus_2);

	hCosThetaK_Y2_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExEt_Y2_p8Be0_3->SetDirectory(p8Be0plus_3);
	hEyEt_Y2_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExEt_CosThetaK_Y2_p8Be0_3->SetDirectory(p8Be0plus_3);
	hEyEt_CosThetaK_Y2_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExEt_Et_Y2_p8Be0_3->SetDirectory(p8Be0plus_3);
	hExSPS_ExEt_Y2_p8Be0_3->SetDirectory(p8Be0plus_3);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//p+8Be(2+)////////////////////////////////////////////////////////////////////////////////////////////////////
	TH2D *hDalitzInvMass_p8Be2_1 = new TH2D("hDalitzInvMass_p8Be2_1", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);
	TH2D *hDalitzInvMass_p8Be2_2 = new TH2D("hDalitzInvMass_p8Be2_2", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);
	TH2D *hDalitzInvMass_p8Be2_3 = new TH2D("hDalitzInvMass_p8Be2_3", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);


	//SABREsum histos
	TH2D *hSABREsumE_vs_ExSPS_p8Be2_1 = new TH2D("hSABREsumE_vs_ExSPS_p8Be2_1","hSABREsumE_vs_ExSPS_p8Be2_1", 1400, 0, 7, 1400, 0, 7);
	TH2D *hSABREsumE_vs_ExSPS_p8Be2_2 = new TH2D("hSABREsumE_vs_ExSPS_p8Be2_2","hSABREsumE_vs_ExSPS_p8Be2_2", 1400, 0, 7, 1400, 0, 7);
	TH2D *hSABREsumE_vs_ExSPS_p8Be2_3 = new TH2D("hSABREsumE_vs_ExSPS_p8Be2_3","hSABREsumE_vs_ExSPS_p8Be2_3", 1400, 0, 7, 1400, 0, 7);

	//9B Ex histos:
	TH1D *hExSPS_p8Be2_1 = new TH1D("hExSPS_p8Be2_1", "hExSPS_p8Be2_1", 1600, -1, 7);
	TH1D *hRecEx_p8Be2_1 = new TH1D("hRecEx_p8Be2_1", "hRecEx_p8Be2_1", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_p8Be2_1 = new TH2D("hRecEx_VS_ExSPS_p8Be2_1", "hRecEx_VS_ExSPS_p8Be2_1", 1600, -1, 7, 1600, -1, 7);

	TH1D *hExSPS_p8Be2_2 = new TH1D("hExSPS_p8Be2_2", "hExSPS_p8Be2_2", 1600, -1, 7);
	TH1D *hRecEx_p8Be2_2 = new TH1D("hRecEx_p8Be2_2", "hRecEx_p8Be2_2", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_p8Be2_2 = new TH2D("hRecEx_VS_ExSPS_p8Be2_2", "hRecEx_VS_ExSPS_p8Be2_2", 1600, -1, 7, 1600, -1, 7);

	TH1D *hExSPS_p8Be2_3 = new TH1D("hExSPS_p8Be2_3", "hExSPS_p8Be2_3", 1600, -1, 7);
	TH1D *hRecEx_p8Be2_3 = new TH1D("hRecEx_p8Be2_3", "hRecEx_p8Be2_3", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_p8Be2_3 = new TH2D("hRecEx_VS_ExSPS_p8Be2_3", "hRecEx_VS_ExSPS_p8Be2_3", 1600, -1, 7, 1600, -1, 7);

	//jacobi histos
		//1sig
			//T
				TH1D *hCosThetaK_T_p8Be2_1 = new TH1D("hCosThetaK_T_p8Be2_1","hCosThetaK_T_p8Be2_1", 200, -1, 1);
				TH1D *hExEt_T_p8Be2_1 = new TH1D("hExEt_T_p8Be2_1","hExEt_T_p8Be2_1",100, 0, 1);
				TH1D *hEyEt_T_p8Be2_1 = new TH1D("hEyEt_T_p8Be2_1", "hEyEt_T_p8Be2_1", 100, 0, 1);
				TH2D *hExEt_CosThetaK_T_p8Be2_1 = new TH2D("hExEt_CosThetaK_T_p8Be2_1","hExEt_CosThetaK_T_p8Be2_1",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_T_p8Be2_1 = new TH2D("hEyEt_CosThetaK_T_p8Be2_1","hEyEt_CosThetaK_T_p8Be2_1",200,-1,1,100,0,1);
				TH2D *hExEt_Et_T_p8Be2_1 = new TH2D("hExEt_Et_T_p8Be2_1", "hExEt_Et_T_p8Be2_1", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_T_p8Be2_1 = new TH2D("hExSPS_ExEt_T_p8Be2_1","hExSPS_ExEt_T_p8Be2_1",100,0,1,1400,0,7);
			//Y1
				TH1D *hCosThetaK_Y1_p8Be2_1 = new TH1D("hCosThetaK_Y1_p8Be2_1","hCosThetaK_Y1_p8Be2_1", 200, -1, 1);
				TH1D *hExEt_Y1_p8Be2_1 = new TH1D("hExEt_Y1_p8Be2_1","hExEt_Y1_p8Be2_1",100, 0, 1);
				TH1D *hEyEt_Y1_p8Be2_1 = new TH1D("hEyEt_Y1_p8Be2_1", "hEyEt_Y1_p8Be2_1", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y1_p8Be2_1 = new TH2D("hExEt_CosThetaK_Y1_p8Be2_1","hExEt_CosThetaK_Y1_p8Be2_1",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y1_p8Be2_1 = new TH2D("hEyEt_CosThetaK_Y1_p8Be2_1","hEyEt_CosThetaK_Y1_p8Be2_1",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y1_p8Be2_1 = new TH2D("hExEt_Et_Y1_p8Be2_1", "hExEt_Et_Y1_p8Be2_1", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y1_p8Be2_1 = new TH2D("hExSPS_ExEt_Y1_p8Be2_1","hExSPS_ExEt_Y1_p8Be2_1",100,0,1,1400,0,7);
			//Y2
				TH1D *hCosThetaK_Y2_p8Be2_1 = new TH1D("hCosThetaK_Y2_p8Be2_1","hCosThetaK_Y2_p8Be2_1", 200, -1, 1);
				TH1D *hExEt_Y2_p8Be2_1 = new TH1D("hExEt_Y2_p8Be2_1","hExEt_Y2_p8Be2_1",100, 0, 1);
				TH1D *hEyEt_Y2_p8Be2_1 = new TH1D("hEyEt_Y2_p8Be2_1", "hEyEt_Y2_p8Be2_1", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y2_p8Be2_1 = new TH2D("hExEt_CosThetaK_Y2_p8Be2_1","hExEt_CosThetaK_Y2_p8Be2_1",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y2_p8Be2_1 = new TH2D("hEyEt_CosThetaK_Y2_p8Be2_1","hEyEt_CosThetaK_Y2_p8Be2_1",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y2_p8Be2_1 = new TH2D("hExEt_Et_Y2_p8Be2_1", "hExEt_Et_Y2_p8Be2_1", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y2_p8Be2_1 = new TH2D("hExSPS_ExEt_Y2_p8Be2_1","hExSPS_ExEt_Y2_p8Be2_1",100,0,1,1400,0,7);
		//2sig
			//T
				TH1D *hCosThetaK_T_p8Be2_2 = new TH1D("hCosThetaK_T_p8Be2_2","hCosThetaK_T_p8Be2_2", 200, -1, 1);
				TH1D *hExEt_T_p8Be2_2 = new TH1D("hExEt_T_p8Be2_2","hExEt_T_p8Be2_2",100, 0, 1);
				TH1D *hEyEt_T_p8Be2_2 = new TH1D("hEyEt_T_p8Be2_2", "hEyEt_T_p8Be2_2", 100, 0, 1);
				TH2D *hExEt_CosThetaK_T_p8Be2_2 = new TH2D("hExEt_CosThetaK_T_p8Be2_2","hExEt_CosThetaK_T_p8Be2_2",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_T_p8Be2_2 = new TH2D("hEyEt_CosThetaK_T_p8Be2_2","hEyEt_CosThetaK_T_p8Be2_2",200,-1,1,100,0,1);
				TH2D *hExEt_Et_T_p8Be2_2 = new TH2D("hExEt_Et_T_p8Be2_2", "hExEt_Et_T_p8Be2_2", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_T_p8Be2_2 = new TH2D("hExSPS_ExEt_T_p8Be2_2","hExSPS_ExEt_T_p8Be2_2",100,0,1,1400,0,7);
			//Y1
				TH1D *hCosThetaK_Y1_p8Be2_2 = new TH1D("hCosThetaK_Y1_p8Be2_2","hCosThetaK_Y1_p8Be2_2", 200, -1, 1);
				TH1D *hExEt_Y1_p8Be2_2 = new TH1D("hExEt_Y1_p8Be2_2","hExEt_Y1_p8Be2_2",100, 0, 1);
				TH1D *hEyEt_Y1_p8Be2_2 = new TH1D("hEyEt_Y1_p8Be2_2", "hEyEt_Y1_p8Be2_2", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y1_p8Be2_2 = new TH2D("hExEt_CosThetaK_Y1_p8Be2_2","hExEt_CosThetaK_Y1_p8Be2_2",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y1_p8Be2_2 = new TH2D("hEyEt_CosThetaK_Y1_p8Be2_2","hEyEt_CosThetaK_Y1_p8Be2_2",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y1_p8Be2_2 = new TH2D("hExEt_Et_Y1_p8Be2_2", "hExEt_Et_Y1_p8Be2_2", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y1_p8Be2_2 = new TH2D("hExSPS_ExEt_Y1_p8Be2_2","hExSPS_ExEt_Y1_p8Be2_2",100,0,1,1400,0,7);
			//Y2
				TH1D *hCosThetaK_Y2_p8Be2_2 = new TH1D("hCosThetaK_Y2_p8Be2_2","hCosThetaK_Y2_p8Be2_2", 200, -1, 1);
				TH1D *hExEt_Y2_p8Be2_2 = new TH1D("hExEt_Y2_p8Be2_2","hExEt_Y2_p8Be2_2",100, 0, 1);
				TH1D *hEyEt_Y2_p8Be2_2 = new TH1D("hEyEt_Y2_p8Be2_2", "hEyEt_Y2_p8Be2_2", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y2_p8Be2_2 = new TH2D("hExEt_CosThetaK_Y2_p8Be2_2","hExEt_CosThetaK_Y2_p8Be2_2",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y2_p8Be2_2 = new TH2D("hEyEt_CosThetaK_Y2_p8Be2_2","hEyEt_CosThetaK_Y2_p8Be2_2",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y2_p8Be2_2 = new TH2D("hExEt_Et_Y2_p8Be2_2", "hExEt_Et_Y2_p8Be2_2", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y2_p8Be2_2 = new TH2D("hExSPS_ExEt_Y2_p8Be2_2","hExSPS_ExEt_Y2_p8Be2_2",100,0,1,1400,0,7);
		//3sig
			//T
				TH1D *hCosThetaK_T_p8Be2_3 = new TH1D("hCosThetaK_T_p8Be2_3","hCosThetaK_T_p8Be2_3", 200, -1, 1);
				TH1D *hExEt_T_p8Be2_3 = new TH1D("hExEt_T_p8Be2_3","hExEt_T_p8Be2_3",100, 0, 1);
				TH1D *hEyEt_T_p8Be2_3 = new TH1D("hEyEt_T_p8Be2_3", "hEyEt_T_p8Be2_3", 100, 0, 1);
				TH2D *hExEt_CosThetaK_T_p8Be2_3 = new TH2D("hExEt_CosThetaK_T_p8Be2_3","hExEt_CosThetaK_T_p8Be2_3",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_T_p8Be2_3 = new TH2D("hEyEt_CosThetaK_T_p8Be2_3","hEyEt_CosThetaK_T_p8Be2_3",200,-1,1,100,0,1);
				TH2D *hExEt_Et_T_p8Be2_3 = new TH2D("hExEt_Et_T_p8Be2_3", "hExEt_Et_T_p8Be2_3", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_T_p8Be2_3 = new TH2D("hExSPS_ExEt_T_p8Be2_3","hExSPS_ExEt_T_p8Be2_3",100,0,1,1400,0,7);
			//Y1
				TH1D *hCosThetaK_Y1_p8Be2_3 = new TH1D("hCosThetaK_Y1_p8Be2_3","hCosThetaK_Y1_p8Be2_3", 200, -1, 1);
				TH1D *hExEt_Y1_p8Be2_3 = new TH1D("hExEt_Y1_p8Be2_3","hExEt_Y1_p8Be2_3",100, 0, 1);
				TH1D *hEyEt_Y1_p8Be2_3 = new TH1D("hEyEt_Y1_p8Be2_3", "hEyEt_Y1_p8Be2_3", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y1_p8Be2_3 = new TH2D("hExEt_CosThetaK_Y1_p8Be2_3","hExEt_CosThetaK_Y1_p8Be2_3",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y1_p8Be2_3 = new TH2D("hEyEt_CosThetaK_Y1_p8Be2_3","hEyEt_CosThetaK_Y1_p8Be2_3",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y1_p8Be2_3 = new TH2D("hExEt_Et_Y1_p8Be2_3", "hExEt_Et_Y1_p8Be2_3", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y1_p8Be2_3 = new TH2D("hExSPS_ExEt_Y1_p8Be2_3","hExSPS_ExEt_Y1_p8Be2_3",100,0,1,1400,0,7);
			//Y2
				TH1D *hCosThetaK_Y2_p8Be2_3 = new TH1D("hCosThetaK_Y2_p8Be2_3","hCosThetaK_Y2_p8Be2_3", 200, -1, 1);
				TH1D *hExEt_Y2_p8Be2_3 = new TH1D("hExEt_Y2_p8Be2_3","hExEt_Y2_p8Be2_3",100, 0, 1);
				TH1D *hEyEt_Y2_p8Be2_3 = new TH1D("hEyEt_Y2_p8Be2_3", "hEyEt_Y2_p8Be2_3", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y2_p8Be2_3 = new TH2D("hExEt_CosThetaK_Y2_p8Be2_3","hExEt_CosThetaK_Y2_p8Be2_3",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y2_p8Be2_3 = new TH2D("hEyEt_CosThetaK_Y2_p8Be2_3","hEyEt_CosThetaK_Y2_p8Be2_3",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y2_p8Be2_3 = new TH2D("hExEt_Et_Y2_p8Be2_3", "hExEt_Et_Y2_p8Be2_3", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y2_p8Be2_3 = new TH2D("hExSPS_ExEt_Y2_p8Be2_3","hExSPS_ExEt_Y2_p8Be2_3",100,0,1,1400,0,7);

	//prep directories:
	TDirectory *p8Be2plus = outfile->mkdir("p8Be2plus");
	TDirectory *p8Be2plus_1 = p8Be2plus->mkdir("p8Be2_1sig");
	TDirectory *p8Be2plus_2 = p8Be2plus->mkdir("p8Be2_2sig");
	TDirectory *p8Be2plus_3 = p8Be2plus->mkdir("p8Be2_3sig");

	hDalitzInvMass_p8Be2_1->SetDirectory(p8Be2plus_1);
	hDalitzInvMass_p8Be2_2->SetDirectory(p8Be2plus_2);
	hDalitzInvMass_p8Be2_3->SetDirectory(p8Be2plus_3);

	hSABREsumE_vs_ExSPS_p8Be2_1->SetDirectory(p8Be2plus_1);
	hSABREsumE_vs_ExSPS_p8Be2_2->SetDirectory(p8Be2plus_2);
	hSABREsumE_vs_ExSPS_p8Be2_3->SetDirectory(p8Be2plus_3);

	hExSPS_p8Be2_1->SetDirectory(p8Be2plus_1);
	hRecEx_p8Be2_1->SetDirectory(p8Be2plus_1);
	hRecEx_VS_ExSPS_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExSPS_p8Be2_2->SetDirectory(p8Be2plus_2);
	hRecEx_p8Be2_2->SetDirectory(p8Be2plus_2);
	hRecEx_VS_ExSPS_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExSPS_p8Be2_3->SetDirectory(p8Be2plus_3);
	hRecEx_p8Be2_3->SetDirectory(p8Be2plus_3);
	hRecEx_VS_ExSPS_p8Be2_3->SetDirectory(p8Be2plus_3);

	hCosThetaK_T_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExEt_T_p8Be2_1->SetDirectory(p8Be2plus_1);
	hEyEt_T_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExEt_CosThetaK_T_p8Be2_1->SetDirectory(p8Be2plus_1);
	hEyEt_CosThetaK_T_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExEt_Et_T_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExSPS_ExEt_T_p8Be2_1->SetDirectory(p8Be2plus_1);

	hCosThetaK_T_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExEt_T_p8Be2_2->SetDirectory(p8Be2plus_2);
	hEyEt_T_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExEt_CosThetaK_T_p8Be2_2->SetDirectory(p8Be2plus_2);
	hEyEt_CosThetaK_T_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExEt_Et_T_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExSPS_ExEt_T_p8Be2_2->SetDirectory(p8Be2plus_2);

	hCosThetaK_T_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExEt_T_p8Be2_3->SetDirectory(p8Be2plus_3);
	hEyEt_T_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExEt_CosThetaK_T_p8Be2_3->SetDirectory(p8Be2plus_3);
	hEyEt_CosThetaK_T_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExEt_Et_T_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExSPS_ExEt_T_p8Be2_3->SetDirectory(p8Be2plus_3);

	hCosThetaK_Y1_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExEt_Y1_p8Be2_1->SetDirectory(p8Be2plus_1);
	hEyEt_Y1_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExEt_CosThetaK_Y1_p8Be2_1->SetDirectory(p8Be2plus_1);
	hEyEt_CosThetaK_Y1_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExEt_Et_Y1_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExSPS_ExEt_Y1_p8Be2_1->SetDirectory(p8Be2plus_1);

	hCosThetaK_Y1_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExEt_Y1_p8Be2_2->SetDirectory(p8Be2plus_2);
	hEyEt_Y1_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExEt_CosThetaK_Y1_p8Be2_2->SetDirectory(p8Be2plus_2);
	hEyEt_CosThetaK_Y1_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExEt_Et_Y1_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExSPS_ExEt_Y1_p8Be2_2->SetDirectory(p8Be2plus_2);

	hCosThetaK_Y1_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExEt_Y1_p8Be2_3->SetDirectory(p8Be2plus_3);
	hEyEt_Y1_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExEt_CosThetaK_Y1_p8Be2_3->SetDirectory(p8Be2plus_3);
	hEyEt_CosThetaK_Y1_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExEt_Et_Y1_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExSPS_ExEt_Y1_p8Be2_3->SetDirectory(p8Be2plus_3);

	hCosThetaK_Y2_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExEt_Y2_p8Be2_1->SetDirectory(p8Be2plus_1);
	hEyEt_Y2_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExEt_CosThetaK_Y2_p8Be2_1->SetDirectory(p8Be2plus_1);
	hEyEt_CosThetaK_Y2_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExEt_Et_Y2_p8Be2_1->SetDirectory(p8Be2plus_1);
	hExSPS_ExEt_Y2_p8Be2_1->SetDirectory(p8Be2plus_1);

	hCosThetaK_Y2_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExEt_Y2_p8Be2_2->SetDirectory(p8Be2plus_2);
	hEyEt_Y2_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExEt_CosThetaK_Y2_p8Be2_2->SetDirectory(p8Be2plus_2);
	hEyEt_CosThetaK_Y2_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExEt_Et_Y2_p8Be2_2->SetDirectory(p8Be2plus_2);
	hExSPS_ExEt_Y2_p8Be2_2->SetDirectory(p8Be2plus_2);

	hCosThetaK_Y2_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExEt_Y2_p8Be2_3->SetDirectory(p8Be2plus_3);
	hEyEt_Y2_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExEt_CosThetaK_Y2_p8Be2_3->SetDirectory(p8Be2plus_3);
	hEyEt_CosThetaK_Y2_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExEt_Et_Y2_p8Be2_3->SetDirectory(p8Be2plus_3);
	hExSPS_ExEt_Y2_p8Be2_3->SetDirectory(p8Be2plus_3);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//a+5Li(gs+)////////////////////////////////////////////////////////////////////////////////////////////////////
	TH2D *hDalitzInvMass_a5Ligs_1 = new TH2D("hDalitzInvMass_a5Ligs_1", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);
	TH2D *hDalitzInvMass_a5Ligs_2 = new TH2D("hDalitzInvMass_a5Ligs_2", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);
	TH2D *hDalitzInvMass_a5Ligs_3 = new TH2D("hDalitzInvMass_a5Ligs_3", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);

	//SABREsum histos
	TH2D *hSABREsumE_vs_ExSPS_a5Ligs_1 = new TH2D("hSABREsumE_vs_ExSPS_a5Ligs_1","hSABREsumE_vs_ExSPS_a5Ligs_1", 1400, 0, 7, 1400, 0, 7);
	TH2D *hSABREsumE_vs_ExSPS_a5Ligs_2 = new TH2D("hSABREsumE_vs_ExSPS_a5Ligs_2","hSABREsumE_vs_ExSPS_a5Ligs_2", 1400, 0, 7, 1400, 0, 7);
	TH2D *hSABREsumE_vs_ExSPS_a5Ligs_3 = new TH2D("hSABREsumE_vs_ExSPS_a5Ligs_3","hSABREsumE_vs_ExSPS_a5Ligs_3", 1400, 0, 7, 1400, 0, 7);

	//9B Ex histos:
	TH1D *hExSPS_a5Ligs_1 = new TH1D("hExSPS_a5Ligs_1", "hExSPS_a5Ligs_1", 1600, -1, 7);
	TH1D *hRecEx_a5Ligs_1 = new TH1D("hRecEx_a5Ligs_1", "hRecEx_a5Ligs_1", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_a5Ligs_1 = new TH2D("hRecEx_VS_ExSPS_a5Ligs_1", "hRecEx_VS_ExSPS_a5Ligs_1", 1600, -1, 7, 1600, -1, 7);

	TH1D *hExSPS_a5Ligs_2 = new TH1D("hExSPS_a5Ligs_2", "hExSPS_a5Ligs_2", 1600, -1, 7);
	TH1D *hRecEx_a5Ligs_2 = new TH1D("hRecEx_a5Ligs_2", "hRecEx_a5Ligs_2", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_a5Ligs_2 = new TH2D("hRecEx_VS_ExSPS_a5Ligs_2", "hRecEx_VS_ExSPS_a5Ligs_2", 1600, -1, 7, 1600, -1, 7);

	TH1D *hExSPS_a5Ligs_3 = new TH1D("hExSPS_a5Ligs_3", "hExSPS_a5Ligs_3", 1600, -1, 7);
	TH1D *hRecEx_a5Ligs_3 = new TH1D("hRecEx_a5Ligs_3", "hRecEx_a5Ligs_3", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_a5Ligs_3 = new TH2D("hRecEx_VS_ExSPS_a5Ligs_3", "hRecEx_VS_ExSPS_a5Ligs_3", 1600, -1, 7, 1600, -1, 7);

	//jacobi histos
	//1sig
		//T
			TH1D *hCosThetaK_T_a5Ligs_1 = new TH1D("hCosThetaK_T_a5Ligs_1","hCosThetaK_T_a5Ligs_1", 200, -1, 1);
			TH1D *hExEt_T_a5Ligs_1 = new TH1D("hExEt_T_a5Ligs_1","hExEt_T_a5Ligs_1",100, 0, 1);
			TH1D *hEyEt_T_a5Ligs_1 = new TH1D("hEyEt_T_a5Ligs_1", "hEyEt_T_a5Ligs_1", 100, 0, 1);
			TH2D *hExEt_CosThetaK_T_a5Ligs_1 = new TH2D("hExEt_CosThetaK_T_a5Ligs_1","hExEt_CosThetaK_T_a5Ligs_1",200,-1,1,100,0,1);
			TH2D *hEyEt_CosThetaK_T_a5Ligs_1 = new TH2D("hEyEt_CosThetaK_T_a5Ligs_1","hEyEt_CosThetaK_T_a5Ligs_1",200,-1,1,100,0,1);
			TH2D *hExEt_Et_T_a5Ligs_1 = new TH2D("hExEt_Et_T_a5Ligs_1", "hExEt_Et_T_a5Ligs_1", 1000, 0, 10, 100, 0, 1);
			TH2D *hExSPS_ExEt_T_a5Ligs_1 = new TH2D("hExSPS_ExEt_T_a5Ligs_1","hExSPS_ExEt_T_a5Ligs_1",100,0,1,1400,0,7);
		//Y1
			TH1D *hCosThetaK_Y1_a5Ligs_1 = new TH1D("hCosThetaK_Y1_a5Ligs_1","hCosThetaK_Y1_a5Ligs_1", 200, -1, 1);
			TH1D *hExEt_Y1_a5Ligs_1 = new TH1D("hExEt_Y1_a5Ligs_1","hExEt_Y1_a5Ligs_1",100, 0, 1);
			TH1D *hEyEt_Y1_a5Ligs_1 = new TH1D("hEyEt_Y1_a5Ligs_1", "hEyEt_Y1_a5Ligs_1", 100, 0, 1);
			TH2D *hExEt_CosThetaK_Y1_a5Ligs_1 = new TH2D("hExEt_CosThetaK_Y1_a5Ligs_1","hExEt_CosThetaK_Y1_a5Ligs_1",200,-1,1,100,0,1);
			TH2D *hEyEt_CosThetaK_Y1_a5Ligs_1 = new TH2D("hEyEt_CosThetaK_Y1_a5Ligs_1","hEyEt_CosThetaK_Y1_a5Ligs_1",200,-1,1,100,0,1);
			TH2D *hExEt_Et_Y1_a5Ligs_1 = new TH2D("hExEt_Et_Y1_a5Ligs_1", "hExEt_Et_Y1_a5Ligs_1", 1000, 0, 10, 100, 0, 1);
			TH2D *hExSPS_ExEt_Y1_a5Ligs_1 = new TH2D("hExSPS_ExEt_Y1_a5Ligs_1","hExSPS_ExEt_Y1_a5Ligs_1",100,0,1,1400,0,7);
		//Y2
			TH1D *hCosThetaK_Y2_a5Ligs_1 = new TH1D("hCosThetaK_Y2_a5Ligs_1","hCosThetaK_Y2_a5Ligs_1", 200, -1, 1);
			TH1D *hExEt_Y2_a5Ligs_1 = new TH1D("hExEt_Y2_a5Ligs_1","hExEt_Y2_a5Ligs_1",100, 0, 1);
			TH1D *hEyEt_Y2_a5Ligs_1 = new TH1D("hEyEt_Y2_a5Ligs_1", "hEyEt_Y2_a5Ligs_1", 100, 0, 1);
			TH2D *hExEt_CosThetaK_Y2_a5Ligs_1 = new TH2D("hExEt_CosThetaK_Y2_a5Ligs_1","hExEt_CosThetaK_Y2_a5Ligs_1",200,-1,1,100,0,1);
			TH2D *hEyEt_CosThetaK_Y2_a5Ligs_1 = new TH2D("hEyEt_CosThetaK_Y2_a5Ligs_1","hEyEt_CosThetaK_Y2_a5Ligs_1",200,-1,1,100,0,1);
			TH2D *hExEt_Et_Y2_a5Ligs_1 = new TH2D("hExEt_Et_Y2_a5Ligs_1", "hExEt_Et_Y2_a5Ligs_1", 1000, 0, 10, 100, 0, 1);
			TH2D *hExSPS_ExEt_Y2_a5Ligs_1 = new TH2D("hExSPS_ExEt_Y2_a5Ligs_1","hExSPS_ExEt_Y2_a5Ligs_1",100,0,1,1400,0,7);
	//2sig
		//T
			TH1D *hCosThetaK_T_a5Ligs_2 = new TH1D("hCosThetaK_T_a5Ligs_2","hCosThetaK_T_a5Ligs_2", 200, -1, 1);
			TH1D *hExEt_T_a5Ligs_2 = new TH1D("hExEt_T_a5Ligs_2","hExEt_T_a5Ligs_2",100, 0, 1);
			TH1D *hEyEt_T_a5Ligs_2 = new TH1D("hEyEt_T_a5Ligs_2", "hEyEt_T_a5Ligs_2", 100, 0, 1);
			TH2D *hExEt_CosThetaK_T_a5Ligs_2 = new TH2D("hExEt_CosThetaK_T_a5Ligs_2","hExEt_CosThetaK_T_a5Ligs_2",200,-1,1,100,0,1);
			TH2D *hEyEt_CosThetaK_T_a5Ligs_2 = new TH2D("hEyEt_CosThetaK_T_a5Ligs_2","hEyEt_CosThetaK_T_a5Ligs_2",200,-1,1,100,0,1);
			TH2D *hExEt_Et_T_a5Ligs_2 = new TH2D("hExEt_Et_T_a5Ligs_2", "hExEt_Et_T_a5Ligs_2", 1000, 0, 10, 100, 0, 1);
			TH2D *hExSPS_ExEt_T_a5Ligs_2 = new TH2D("hExSPS_ExEt_T_a5Ligs_2","hExSPS_ExEt_T_a5Ligs_2",100,0,1,1400,0,7);
		//Y1
			TH1D *hCosThetaK_Y1_a5Ligs_2 = new TH1D("hCosThetaK_Y1_a5Ligs_2","hCosThetaK_Y1_a5Ligs_2", 200, -1, 1);
			TH1D *hExEt_Y1_a5Ligs_2 = new TH1D("hExEt_Y1_a5Ligs_2","hExEt_Y1_a5Ligs_2",100, 0, 1);
			TH1D *hEyEt_Y1_a5Ligs_2 = new TH1D("hEyEt_Y1_a5Ligs_2", "hEyEt_Y1_a5Ligs_2", 100, 0, 1);
			TH2D *hExEt_CosThetaK_Y1_a5Ligs_2 = new TH2D("hExEt_CosThetaK_Y1_a5Ligs_2","hExEt_CosThetaK_Y1_a5Ligs_2",200,-1,1,100,0,1);
			TH2D *hEyEt_CosThetaK_Y1_a5Ligs_2 = new TH2D("hEyEt_CosThetaK_Y1_a5Ligs_2","hEyEt_CosThetaK_Y1_a5Ligs_2",200,-1,1,100,0,1);
			TH2D *hExEt_Et_Y1_a5Ligs_2 = new TH2D("hExEt_Et_Y1_a5Ligs_2", "hExEt_Et_Y1_a5Ligs_2", 1000, 0, 10, 100, 0, 1);
			TH2D *hExSPS_ExEt_Y1_a5Ligs_2 = new TH2D("hExSPS_ExEt_Y1_a5Ligs_2","hExSPS_ExEt_Y1_a5Ligs_2",100,0,1,1400,0,7);
		//Y2
			TH1D *hCosThetaK_Y2_a5Ligs_2 = new TH1D("hCosThetaK_Y2_a5Ligs_2","hCosThetaK_Y2_a5Ligs_2", 200, -1, 1);
			TH1D *hExEt_Y2_a5Ligs_2 = new TH1D("hExEt_Y2_a5Ligs_2","hExEt_Y2_a5Ligs_2",100, 0, 1);
			TH1D *hEyEt_Y2_a5Ligs_2 = new TH1D("hEyEt_Y2_a5Ligs_2", "hEyEt_Y2_a5Ligs_2", 100, 0, 1);
			TH2D *hExEt_CosThetaK_Y2_a5Ligs_2 = new TH2D("hExEt_CosThetaK_Y2_a5Ligs_2","hExEt_CosThetaK_Y2_a5Ligs_2",200,-1,1,100,0,1);
			TH2D *hEyEt_CosThetaK_Y2_a5Ligs_2 = new TH2D("hEyEt_CosThetaK_Y2_a5Ligs_2","hEyEt_CosThetaK_Y2_a5Ligs_2",200,-1,1,100,0,1);
			TH2D *hExEt_Et_Y2_a5Ligs_2 = new TH2D("hExEt_Et_Y2_a5Ligs_2", "hExEt_Et_Y2_a5Ligs_2", 1000, 0, 10, 100, 0, 1);
			TH2D *hExSPS_ExEt_Y2_a5Ligs_2 = new TH2D("hExSPS_ExEt_Y2_a5Ligs_2","hExSPS_ExEt_Y2_a5Ligs_2",100,0,1,1400,0,7);
	//3sig
		//T
			TH1D *hCosThetaK_T_a5Ligs_3 = new TH1D("hCosThetaK_T_a5Ligs_3","hCosThetaK_T_a5Ligs_3", 200, -1, 1);
			TH1D *hExEt_T_a5Ligs_3 = new TH1D("hExEt_T_a5Ligs_3","hExEt_T_a5Ligs_3",100, 0, 1);
			TH1D *hEyEt_T_a5Ligs_3 = new TH1D("hEyEt_T_a5Ligs_3", "hEyEt_T_a5Ligs_3", 100, 0, 1);
			TH2D *hExEt_CosThetaK_T_a5Ligs_3 = new TH2D("hExEt_CosThetaK_T_a5Ligs_3","hExEt_CosThetaK_T_a5Ligs_3",200,-1,1,100,0,1);
			TH2D *hEyEt_CosThetaK_T_a5Ligs_3 = new TH2D("hEyEt_CosThetaK_T_a5Ligs_3","hEyEt_CosThetaK_T_a5Ligs_3",200,-1,1,100,0,1);
			TH2D *hExEt_Et_T_a5Ligs_3 = new TH2D("hExEt_Et_T_a5Ligs_3", "hExEt_Et_T_a5Ligs_3", 1000, 0, 10, 100, 0, 1);
			TH2D *hExSPS_ExEt_T_a5Ligs_3 = new TH2D("hExSPS_ExEt_T_a5Ligs_3","hExSPS_ExEt_T_a5Ligs_3",100,0,1,1400,0,7);
		//Y1
			TH1D *hCosThetaK_Y1_a5Ligs_3 = new TH1D("hCosThetaK_Y1_a5Ligs_3","hCosThetaK_Y1_a5Ligs_3", 200, -1, 1);
			TH1D *hExEt_Y1_a5Ligs_3 = new TH1D("hExEt_Y1_a5Ligs_3","hExEt_Y1_a5Ligs_3",100, 0, 1);
			TH1D *hEyEt_Y1_a5Ligs_3 = new TH1D("hEyEt_Y1_a5Ligs_3", "hEyEt_Y1_a5Ligs_3", 100, 0, 1);
			TH2D *hExEt_CosThetaK_Y1_a5Ligs_3 = new TH2D("hExEt_CosThetaK_Y1_a5Ligs_3","hExEt_CosThetaK_Y1_a5Ligs_3",200,-1,1,100,0,1);
			TH2D *hEyEt_CosThetaK_Y1_a5Ligs_3 = new TH2D("hEyEt_CosThetaK_Y1_a5Ligs_3","hEyEt_CosThetaK_Y1_a5Ligs_3",200,-1,1,100,0,1);
			TH2D *hExEt_Et_Y1_a5Ligs_3 = new TH2D("hExEt_Et_Y1_a5Ligs_3", "hExEt_Et_Y1_a5Ligs_3", 1000, 0, 10, 100, 0, 1);
			TH2D *hExSPS_ExEt_Y1_a5Ligs_3 = new TH2D("hExSPS_ExEt_Y1_a5Ligs_3","hExSPS_ExEt_Y1_a5Ligs_3",100,0,1,1400,0,7);
		//Y2
			TH1D *hCosThetaK_Y2_a5Ligs_3 = new TH1D("hCosThetaK_Y2_a5Ligs_3","hCosThetaK_Y2_a5Ligs_3", 200, -1, 1);
			TH1D *hExEt_Y2_a5Ligs_3 = new TH1D("hExEt_Y2_a5Ligs_3","hExEt_Y2_a5Ligs_3",100, 0, 1);
			TH1D *hEyEt_Y2_a5Ligs_3 = new TH1D("hEyEt_Y2_a5Ligs_3", "hEyEt_Y2_a5Ligs_3", 100, 0, 1);
			TH2D *hExEt_CosThetaK_Y2_a5Ligs_3 = new TH2D("hExEt_CosThetaK_Y2_a5Ligs_3","hExEt_CosThetaK_Y2_a5Ligs_3",200,-1,1,100,0,1);
			TH2D *hEyEt_CosThetaK_Y2_a5Ligs_3 = new TH2D("hEyEt_CosThetaK_Y2_a5Ligs_3","hEyEt_CosThetaK_Y2_a5Ligs_3",200,-1,1,100,0,1);
			TH2D *hExEt_Et_Y2_a5Ligs_3 = new TH2D("hExEt_Et_Y2_a5Ligs_3", "hExEt_Et_Y2_a5Ligs_3", 1000, 0, 10, 100, 0, 1);
			TH2D *hExSPS_ExEt_Y2_a5Ligs_3 = new TH2D("hExSPS_ExEt_Y2_a5Ligs_3","hExSPS_ExEt_Y2_a5Ligs_3",100,0,1,1400,0,7);

	//prep directories:
	TDirectory *a5Ligs = outfile->mkdir("a5Ligs");
	TDirectory *a5Ligs_1 = a5Ligs->mkdir("a5Ligs_1sig");
	TDirectory *a5Ligs_2 = a5Ligs->mkdir("a5Ligs_2sig");
	TDirectory *a5Ligs_3 = a5Ligs->mkdir("a5Ligs_3sig");

	hDalitzInvMass_a5Ligs_1->SetDirectory(a5Ligs_1);
	hDalitzInvMass_a5Ligs_2->SetDirectory(a5Ligs_2);
	hDalitzInvMass_a5Ligs_3->SetDirectory(a5Ligs_3);

	hSABREsumE_vs_ExSPS_a5Ligs_1->SetDirectory(a5Ligs_1);
	hSABREsumE_vs_ExSPS_a5Ligs_2->SetDirectory(a5Ligs_2);
	hSABREsumE_vs_ExSPS_a5Ligs_3->SetDirectory(a5Ligs_3);

	hExSPS_a5Ligs_1->SetDirectory(a5Ligs_1);
	hRecEx_a5Ligs_1->SetDirectory(a5Ligs_1);
	hRecEx_VS_ExSPS_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExSPS_a5Ligs_2->SetDirectory(a5Ligs_2);
	hRecEx_a5Ligs_2->SetDirectory(a5Ligs_2);
	hRecEx_VS_ExSPS_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExSPS_a5Ligs_3->SetDirectory(a5Ligs_3);
	hRecEx_a5Ligs_3->SetDirectory(a5Ligs_3);
	hRecEx_VS_ExSPS_a5Ligs_3->SetDirectory(a5Ligs_3);

	hCosThetaK_T_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExEt_T_a5Ligs_1->SetDirectory(a5Ligs_1);
	hEyEt_T_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExEt_CosThetaK_T_a5Ligs_1->SetDirectory(a5Ligs_1);
	hEyEt_CosThetaK_T_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExEt_Et_T_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExSPS_ExEt_T_a5Ligs_1->SetDirectory(a5Ligs_1);

	hCosThetaK_T_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExEt_T_a5Ligs_2->SetDirectory(a5Ligs_2);
	hEyEt_T_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExEt_CosThetaK_T_a5Ligs_2->SetDirectory(a5Ligs_2);
	hEyEt_CosThetaK_T_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExEt_Et_T_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExSPS_ExEt_T_a5Ligs_2->SetDirectory(a5Ligs_2);

	hCosThetaK_T_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExEt_T_a5Ligs_3->SetDirectory(a5Ligs_3);
	hEyEt_T_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExEt_CosThetaK_T_a5Ligs_3->SetDirectory(a5Ligs_3);
	hEyEt_CosThetaK_T_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExEt_Et_T_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExSPS_ExEt_T_a5Ligs_3->SetDirectory(a5Ligs_3);

	hCosThetaK_Y1_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExEt_Y1_a5Ligs_1->SetDirectory(a5Ligs_1);
	hEyEt_Y1_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExEt_CosThetaK_Y1_a5Ligs_1->SetDirectory(a5Ligs_1);
	hEyEt_CosThetaK_Y1_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExEt_Et_Y1_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExSPS_ExEt_Y1_a5Ligs_1->SetDirectory(a5Ligs_1);

	hCosThetaK_Y1_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExEt_Y1_a5Ligs_2->SetDirectory(a5Ligs_2);
	hEyEt_Y1_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExEt_CosThetaK_Y1_a5Ligs_2->SetDirectory(a5Ligs_2);
	hEyEt_CosThetaK_Y1_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExEt_Et_Y1_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExSPS_ExEt_Y1_a5Ligs_2->SetDirectory(a5Ligs_2);

	hCosThetaK_Y1_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExEt_Y1_a5Ligs_3->SetDirectory(a5Ligs_3);
	hEyEt_Y1_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExEt_CosThetaK_Y1_a5Ligs_3->SetDirectory(a5Ligs_3);
	hEyEt_CosThetaK_Y1_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExEt_Et_Y1_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExSPS_ExEt_Y1_a5Ligs_3->SetDirectory(a5Ligs_3);

	hCosThetaK_Y2_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExEt_Y2_a5Ligs_1->SetDirectory(a5Ligs_1);
	hEyEt_Y2_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExEt_CosThetaK_Y2_a5Ligs_1->SetDirectory(a5Ligs_1);
	hEyEt_CosThetaK_Y2_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExEt_Et_Y2_a5Ligs_1->SetDirectory(a5Ligs_1);
	hExSPS_ExEt_Y2_a5Ligs_1->SetDirectory(a5Ligs_1);

	hCosThetaK_Y2_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExEt_Y2_a5Ligs_2->SetDirectory(a5Ligs_2);
	hEyEt_Y2_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExEt_CosThetaK_Y2_a5Ligs_2->SetDirectory(a5Ligs_2);
	hEyEt_CosThetaK_Y2_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExEt_Et_Y2_a5Ligs_2->SetDirectory(a5Ligs_2);
	hExSPS_ExEt_Y2_a5Ligs_2->SetDirectory(a5Ligs_2);

	hCosThetaK_Y2_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExEt_Y2_a5Ligs_3->SetDirectory(a5Ligs_3);
	hEyEt_Y2_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExEt_CosThetaK_Y2_a5Ligs_3->SetDirectory(a5Ligs_3);
	hEyEt_CosThetaK_Y2_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExEt_Et_Y2_a5Ligs_3->SetDirectory(a5Ligs_3);
	hExSPS_ExEt_Y2_a5Ligs_3->SetDirectory(a5Ligs_3);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////INCONSISTENT EVENTS //////////////////////////////////////////////////////////////////////////////////////////
	
	//dalitz plots
	TH2D *hDalitzInvMass_none = new TH2D("hDalitzInvMass_none", "M^{2}_{p+#alpha} vs M^{2}_{#alpha+#alpha};M^{2}_{#alpha+#alpha};M^{2}_{p+#alpha}",
											800, 55.57e6, 55.69e6,
											800, 21.76e6, 21.84e6);

	//SABREsum histos
	TH2D *hSABREsumE_vs_ExSPS_none = new TH2D("hSABREsumE_vs_ExSPS_none","hSABREsumE_vs_ExSPS_none", 1400, 0, 7, 1400, 0, 7);

	//9B Ex histos:
	TH1D *hExSPS_none = new TH1D("hExSPS_none", "hExSPS_none", 1600, -1, 7);
	TH1D *hRecEx_none = new TH1D("hRecEx_none", "hRecEx_none", 1600, -1, 7);
	TH2D *hRecEx_VS_ExSPS_none = new TH2D("hRecEx_VS_ExSPS_none", "hRecEx_VS_ExSPS_none", 1600, -1, 7, 1600, -1, 7);

	//jacobi histos
		//1sig
			//T
				TH1D *hCosThetaK_T_none = new TH1D("hCosThetaK_T_none","hCosThetaK_T_none", 200, -1, 1);
				TH1D *hExEt_T_none = new TH1D("hExEt_T_none","hExEt_T_none",100, 0, 1);
				TH1D *hEyEt_T_none = new TH1D("hEyEt_T_none", "hEyEt_T_none", 100, 0, 1);
				TH2D *hExEt_CosThetaK_T_none = new TH2D("hExEt_CosThetaK_T_none","hExEt_CosThetaK_T_none",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_T_none = new TH2D("hEyEt_CosThetaK_T_none","hEyEt_CosThetaK_T_none",200,-1,1,100,0,1);
				TH2D *hExEt_Et_T_none = new TH2D("hExEt_Et_T_none", "hExEt_Et_T_none", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_T_none = new TH2D("hExSPS_ExEt_T_none","hExSPS_ExEt_T_none",100,0,1,1400,0,7);
			//Y1
				TH1D *hCosThetaK_Y1_none = new TH1D("hCosThetaK_Y1_none","hCosThetaK_Y1_none", 200, -1, 1);
				TH1D *hExEt_Y1_none = new TH1D("hExEt_Y1_none","hExEt_Y1_none",100, 0, 1);
				TH1D *hEyEt_Y1_none = new TH1D("hEyEt_Y1_none", "hEyEt_Y1_none", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y1_none = new TH2D("hExEt_CosThetaK_Y1_none","hExEt_CosThetaK_Y1_none",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y1_none = new TH2D("hEyEt_CosThetaK_Y1_none","hEyEt_CosThetaK_Y1_none",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y1_none = new TH2D("hExEt_Et_Y1_none", "hExEt_Et_Y1_none", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y1_none = new TH2D("hExSPS_ExEt_Y1_none","hExSPS_ExEt_Y1_none",100,0,1,1400,0,7);
			//Y2
				TH1D *hCosThetaK_Y2_none = new TH1D("hCosThetaK_Y2_none","hCosThetaK_Y2_none", 200, -1, 1);
				TH1D *hExEt_Y2_none = new TH1D("hExEt_Y2_none","hExEt_Y2_none",100, 0, 1);
				TH1D *hEyEt_Y2_none = new TH1D("hEyEt_Y2_none", "hEyEt_Y2_none", 100, 0, 1);
				TH2D *hExEt_CosThetaK_Y2_none = new TH2D("hExEt_CosThetaK_Y2_none","hExEt_CosThetaK_Y2_none",200,-1,1,100,0,1);
				TH2D *hEyEt_CosThetaK_Y2_none = new TH2D("hEyEt_CosThetaK_Y2_none","hEyEt_CosThetaK_Y2_none",200,-1,1,100,0,1);
				TH2D *hExEt_Et_Y2_none = new TH2D("hExEt_Et_Y2_none", "hExEt_Et_Y2_none", 1000, 0, 10, 100, 0, 1);
				TH2D *hExSPS_ExEt_Y2_none = new TH2D("hExSPS_ExEt_Y2_none","hExSPS_ExEt_Y2_none",100,0,1,1400,0,7);

	//prep directories:
	TDirectory *none = outfile->mkdir("none");

	//set histograms owner to proper directory
	hDalitzInvMass_none->SetDirectory(none);


	hSABREsumE_vs_ExSPS_none->SetDirectory(none);


	hExSPS_none->SetDirectory(none);
	hRecEx_none->SetDirectory(none);
	hRecEx_VS_ExSPS_none->SetDirectory(none);

	hCosThetaK_T_none->SetDirectory(none);
	hExEt_T_none->SetDirectory(none);
	hEyEt_T_none->SetDirectory(none);
	hExEt_CosThetaK_T_none->SetDirectory(none);
	hEyEt_CosThetaK_T_none->SetDirectory(none);
	hExEt_Et_T_none->SetDirectory(none);
	hExSPS_ExEt_T_none->SetDirectory(none);

	hCosThetaK_Y1_none->SetDirectory(none);
	hExEt_Y1_none->SetDirectory(none);
	hEyEt_Y1_none->SetDirectory(none);
	hExEt_CosThetaK_Y1_none->SetDirectory(none);
	hEyEt_CosThetaK_Y1_none->SetDirectory(none);
	hExEt_Et_Y1_none->SetDirectory(none);
	hExSPS_ExEt_Y1_none->SetDirectory(none);

	hCosThetaK_Y2_none->SetDirectory(none);
	hExEt_Y2_none->SetDirectory(none);
	hEyEt_Y2_none->SetDirectory(none);
	hExEt_CosThetaK_Y2_none->SetDirectory(none);
	hEyEt_CosThetaK_Y2_none->SetDirectory(none);
	hExEt_Et_Y2_none->SetDirectory(none);
	hExSPS_ExEt_Y2_none->SetDirectory(none);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//prepare TTreeReader:
	TTreeReader reader(tree);

	//prepare TTreeReaderValues:
	TTreeReaderValue<double> m2_aa(reader, "m2_aa");
	TTreeReaderValue<double> m2_pa1(reader, "m2_pa1");
	TTreeReaderValue<double> m2_pa2(reader, "m2_pa2");
	TTreeReaderValue<double> ExSPS(reader, "ExSPS");
	TTreeReaderValue<double> RecEx(reader, "RecEx");
	TTreeReaderValue<double> SABREsumE(reader, "SABREsumE");

	TTreeReaderValue<double> E_xT(reader, "E_xT");//energy associated w/ x jacobi coordinate in jacobi T 
	TTreeReaderValue<double> E_yT(reader, "E_yT");//energy associated w/ y jacobi coordinate in jacobi T
	TTreeReaderValue<double> E_Txy_T(reader, "E_Txy_T");//total energy E_T = E_x + E_y in jacobi T
	TTreeReaderValue<double> costhetakT(reader, "costhetakT");//cos(theta_k) in jacobi T

	TTreeReaderValue<double> E_xY1(reader, "E_xY1");//energy associated w/ x jacobi coordinate in jacobi Y1
	TTreeReaderValue<double> E_yY1(reader, "E_yY1");//energy associated w/ y jacobi coordinate in jacobi Y1
	TTreeReaderValue<double> E_Txy_Y1(reader, "E_Txy_Y1");//total energy E_T = E_x + E_y in jacobi Y1
	TTreeReaderValue<double> costhetakY1(reader, "costhetakY1");//cos(theta_k) in jacobi Y1

	TTreeReaderValue<double> E_xY2(reader, "E_xY2");//energy associated w/ x jacobi coordinate in jacobi Y2
	TTreeReaderValue<double> E_yY2(reader, "E_yY2");//energy associated w/ y jacobi coordinate in jacobi Y2
	TTreeReaderValue<double> E_Txy_Y2(reader, "E_Txy_Y2");//total energy E_T = E_x + E_y in jacobi Y2
	TTreeReaderValue<double> costhetakY2(reader, "costhetakY2");//cos(theta_k) in jacobi Y2

	TTreeReaderValue<double> E_TIM(reader, "E_TIM");//total energy calculated via IM (can compare with E_Txy_# )


	//iterate through TTree entries:
	//std::cout << "about to enter loop" << std::endl;
	while(reader.Next()){
		//std::cout << "entry = " << reader.GetCurrentEntry() << std::endl;
		// std::cout << "m2_aa status: " << m2_aa.GetSetupStatus() << std::endl;
		// std::cout << "m2_pa1 status: " << m2_pa1.GetSetupStatus() << std::endl;
		// std::cout << "m2_pa2 status: " << m2_pa2.GetSetupStatus() << std::endl;
		// std::cout << "ExSPS status: " << ExSPS.GetSetupStatus() << std::endl;
		// std::cout << "RecEx status: " << RecEx.GetSetupStatus() << std::endl;
		// std::cout << "SABREsumE status: " << SABREsumE.GetSetupStatus() << std::endl;
		// std::cout << "E_xT status: " << E_xT.GetSetupStatus() << std::endl;
		// std::cout << "E_yT status: " << E_yT.GetSetupStatus() << std::endl;
		// std::cout << "E_Txy_T status: " << E_Txy_T.GetSetupStatus() << std::endl;
		// std::cout << "costhetakT status: " << costhetakT.GetSetupStatus() << std::endl;
		// std::cout << "E_xY1 status: " << E_xY1.GetSetupStatus() << std::endl;
		// std::cout << "E_yY1 status: " << E_yY1.GetSetupStatus() << std::endl;
		// std::cout << "E_Txy_Y1 status: " << E_Txy_Y1.GetSetupStatus() << std::endl;
		// std::cout << "costhetakY1 status: " << costhetakY1.GetSetupStatus() << std::endl;
		// std::cout << "E_xY2 status: " << E_xY2.GetSetupStatus() << std::endl;
		// std::cout << "E_yY2 status: " << E_yY2.GetSetupStatus() << std::endl;
		// std::cout << "E_Txy_Y2 status: " << E_Txy_Y2.GetSetupStatus() << std::endl;
		// std::cout << "costhetakY2 status: " << costhetakY2.GetSetupStatus() << std::endl;
		// std::cout << "E_TIM status: " << E_TIM.GetSetupStatus() << std::endl;
		// std::cout << "TEST 1: reading m2_aa" << std::endl;
		// double aa = *m2_aa;

		// std::cout << "TEST 2: reading m2_pa1" << std::endl;
		// double pa1 = *m2_pa1;

		// std::cout << "TEST 3: reading m2_pa2" << std::endl;
		// double pa2 = *m2_pa2;

		// std::cout << "TEST 4: histogram pointer" << std::endl;
		// std::cout << (void*)hDalitzInvMass << std::endl;

		// std::cout << "TEST 5: first histogram fill" << std::endl;
		// hDalitzInvMass->Fill(aa, pa1);

		// std::cout << "TEST 6: second histogram fill" << std::endl;
		// hDalitzInvMass->Fill(aa, pa2);

		// std::cout << "TEST 7: completed" << std::endl;
		//first, fill the cumulative histograms:
		hDalitzInvMass->Fill(*m2_aa, *m2_pa1);
		hDalitzInvMass->Fill(*m2_aa, *m2_pa2);
		//std::cout << "dalitz filled" << std::endl;

		hSABREsumE_vs_ExSPS->Fill(*ExSPS, *SABREsumE);
		//std::cout << "sabresum filled" << std::endl;

		hExSPS->Fill(*ExSPS);
		hRecEx->Fill(*RecEx);
		hRecEx_VS_ExSPS->Fill(*ExSPS, *RecEx);
	
		hCosThetaK_T->Fill(*costhetakT);
		hExEt_T->Fill((*E_xT)/(*E_xT + *E_yT));
		hEyEt_T->Fill((*E_yT)/(*E_xT + *E_yT));
		hExEt_CosThetaK_T->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
		hEyEt_CosThetaK_T->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
		hExEt_Et_T->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
		hExSPS_ExEt_T->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);
		//std::cout << "jacobi T filled" << std::endl;

		hCosThetaK_Y1->Fill(*costhetakY1);
		hExEt_Y1->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
		hEyEt_Y1->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
		hExEt_CosThetaK_Y1->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
		hEyEt_CosThetaK_Y1->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
		hExEt_Et_Y1->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
		hExSPS_ExEt_Y1->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);
		//std::cout << "jacobi Y1 filled" << std::endl;

		hCosThetaK_Y2->Fill(*costhetakY2);
		hExEt_Y2->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
		hEyEt_Y2->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
		hExEt_CosThetaK_Y2->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
		hEyEt_CosThetaK_Y2->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
		hExEt_Et_Y2->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
		hExSPS_ExEt_Y2->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		//std::cout << "jacobi Y2 filled" << std::endl;


		//check each gate here:
		p8Be0plus_sig1_check = IsInside(*m2_aa, p8Be0plus_sig1);
		p8Be0plus_sig2_check = IsInside(*m2_aa, p8Be0plus_sig2);
		p8Be0plus_sig3_check = IsInside(*m2_aa, p8Be0plus_sig3);
		//std::cout << "p+8Be(0+) check (1/2/3): " << p8Be0plus_sig1_check << ", " << p8Be0plus_sig2_check << ", " << p8Be0plus_sig3_check << std::endl;

		p8Be2plus_sig1_check = IsInside(*m2_aa, p8Be2plus_sig1);
		p8Be2plus_sig2_check = IsInside(*m2_aa, p8Be2plus_sig2);
		p8Be2plus_sig3_check = IsInside(*m2_aa, p8Be2plus_sig3);
		//std::cout << "p+8Be(2+) check (1/2/3): " << p8Be2plus_sig1_check << ", " << p8Be2plus_sig2_check << ", " << p8Be2plus_sig3_check << std::endl;

		a5Ligs_sig1_check = IsInside(*m2_pa1, a5Ligs_sig1) || IsInside(*m2_pa2, a5Ligs_sig1);
		a5Ligs_sig2_check = IsInside(*m2_pa1, a5Ligs_sig2) || IsInside(*m2_pa2, a5Ligs_sig2);
		a5Ligs_sig3_check = IsInside(*m2_pa1, a5Ligs_sig3) || IsInside(*m2_pa2, a5Ligs_sig3);
		//std::cout << "a+5Li(gs) check (1/2/3): " << a5Ligs_sig1_check << ", " << a5Ligs_sig2_check << ", " << a5Ligs_sig3_check << std::endl;

		int p8Be0plus_sigmaNum = -1, p8Be2plus_sigmaNum = -1, a5Ligs_sigmaNum = -1;

		if(p8Be0plus_sig1_check){
			p8Be0plus_sigmaNum = 1;
		} else if(p8Be0plus_sig2_check){
			p8Be0plus_sigmaNum = 2;
		} else if(p8Be0plus_sig3_check){
			p8Be0plus_sigmaNum = 3;
		}

		if(p8Be2plus_sig1_check){
			p8Be2plus_sigmaNum = 1;
		} else if(p8Be2plus_sig2_check){
			p8Be2plus_sigmaNum = 2;
		} else if(p8Be2plus_sig3_check){
			p8Be2plus_sigmaNum = 3;
		}

		if(a5Ligs_sig1_check){
			a5Ligs_sigmaNum = 1;
		} else if(p8Be2plus_sig2_check){
			a5Ligs_sigmaNum = 2;
		} else if(p8Be2plus_sig3_check){
			a5Ligs_sigmaNum = 3;
		}


		if(a5Ligs_sigmaNum != -1){
			//fill h2_a5Ligs_VS_p8Be0plus
			if(p8Be0plus_sigmaNum != -1) h2_a5Ligs_VS_p8Be0plus->Fill(p8Be0plus_sigmaNum, a5Ligs_sigmaNum);

			//fill h2_a5Ligs_VS_p8Be2plus
			if(p8Be2plus_sigmaNum != -1) h2_a5Ligs_VS_p8Be2plus->Fill(p8Be2plus_sigmaNum, a5Ligs_sigmaNum);
		}


		//next, check p+8Be(0+):
		if(p8Be0plus_sig1_check){
			//std::cout << "entering p8Be0plus_sig1_check" << std::endl;
			//events that are +/- 1 sigma from expected p+8Be(0+) IM2
			//first, fill the cumulative histograms:
			hDalitzInvMass_p8Be0_1->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_p8Be0_1->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_p8Be0_1->Fill(*ExSPS, *SABREsumE);

			hExSPS_p8Be0_1->Fill(*ExSPS);
			hRecEx_p8Be0_1->Fill(*RecEx);
			hRecEx_VS_ExSPS_p8Be0_1->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_p8Be0_1->Fill(*costhetakT);
			hExEt_T_p8Be0_1->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_p8Be0_1->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_p8Be0_1->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_p8Be0_1->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_p8Be0_1->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_p8Be0_1->Fill(*ExSPS, (*E_xT)/(*E_xT + *E_yT));

			hCosThetaK_Y1_p8Be0_1->Fill(*costhetakY1);
			hExEt_Y1_p8Be0_1->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_p8Be0_1->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_p8Be0_1->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_p8Be0_1->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_p8Be0_1->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_p8Be0_1->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_p8Be0_1->Fill(*costhetakY2);
			hExEt_Y2_p8Be0_1->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_p8Be0_1->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_p8Be0_1->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_p8Be0_1->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_p8Be0_1->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_p8Be0_1->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}
		if(p8Be0plus_sig2_check){ //else if would act as exclusive gates, but we want 2sig to include 1sig and 3sig to include 1sig and 2sig
			//events that are +/-2 sigma from expected p+8Be(0+) IM2
			//std::cout << "entering p8Be0plus_sig2_check" << std::endl;
			//first, fill the cumulative histograms:
			hDalitzInvMass_p8Be0_2->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_p8Be0_2->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_p8Be0_2->Fill(*ExSPS, *SABREsumE);

			hExSPS_p8Be0_2->Fill(*ExSPS);
			hRecEx_p8Be0_2->Fill(*RecEx);
			hRecEx_VS_ExSPS_p8Be0_2->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_p8Be0_2->Fill(*costhetakT);
			hExEt_T_p8Be0_2->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_p8Be0_2->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_p8Be0_2->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_p8Be0_2->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_p8Be0_2->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_p8Be0_2->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);

			hCosThetaK_Y1_p8Be0_2->Fill(*costhetakY1);
			hExEt_Y1_p8Be0_2->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_p8Be0_2->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_p8Be0_2->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_p8Be0_2->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_p8Be0_2->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_p8Be0_2->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_p8Be0_2->Fill(*costhetakY2);
			hExEt_Y2_p8Be0_2->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_p8Be0_2->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_p8Be0_2->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_p8Be0_2->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_p8Be0_2->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_p8Be0_2->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}
		if(p8Be0plus_sig3_check){ //else if would act as exclusive gates, but we want 2sig to include 1sig and 3sig to include 1sig and 2sig
			//events that are +/-3 sigma from expected p+8Be(0+) IM2
			//std::cout << "entering p8Be0plus_sig3_check" << std::endl;
			//first, fill the cumulative histograms:
			hDalitzInvMass_p8Be0_3->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_p8Be0_3->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_p8Be0_3->Fill(*ExSPS, *SABREsumE);

			hExSPS_p8Be0_3->Fill(*ExSPS);
			hRecEx_p8Be0_3->Fill(*RecEx);
			hRecEx_VS_ExSPS_p8Be0_3->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_p8Be0_3->Fill(*costhetakT);
			hExEt_T_p8Be0_3->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_p8Be0_3->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_p8Be0_3->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_p8Be0_3->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_p8Be0_3->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_p8Be0_3->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);

			hCosThetaK_Y1_p8Be0_3->Fill(*costhetakY1);
			hExEt_Y1_p8Be0_3->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_p8Be0_3->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_p8Be0_3->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_p8Be0_3->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_p8Be0_3->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_p8Be0_3->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_p8Be0_3->Fill(*costhetakY2);
			hExEt_Y2_p8Be0_3->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_p8Be0_3->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_p8Be0_3->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_p8Be0_3->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_p8Be0_3->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_p8Be0_3->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}

		//now, check p+8Be(2+):
		if(p8Be2plus_sig1_check){
			//events that are +/- 1 sigma from expected p+8Be(2+) IM2
			//std::cout << "entering p8Be2plus_sig1_check" << std::endl;
			//first, fill the cumulative histograms:
			hDalitzInvMass_p8Be2_1->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_p8Be2_1->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_p8Be2_1->Fill(*ExSPS, *SABREsumE);

			hExSPS_p8Be2_1->Fill(*ExSPS);
			hRecEx_p8Be2_1->Fill(*RecEx);
			hRecEx_VS_ExSPS_p8Be2_1->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_p8Be2_1->Fill(*costhetakT);
			hExEt_T_p8Be2_1->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_p8Be2_1->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_p8Be2_1->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_p8Be2_1->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_p8Be2_1->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_p8Be2_1->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);

			hCosThetaK_Y1_p8Be2_1->Fill(*costhetakY1);
			hExEt_Y1_p8Be2_1->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_p8Be2_1->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_p8Be2_1->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_p8Be2_1->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_p8Be2_1->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_p8Be2_1->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_p8Be2_1->Fill(*costhetakY2);
			hExEt_Y2_p8Be2_1->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_p8Be2_1->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_p8Be2_1->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_p8Be2_1->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_p8Be2_1->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_p8Be2_1->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}
		if(p8Be2plus_sig2_check){ //else if would act as exclusive gates, but we want 2sig to include 1sig and 3sig to include 1sig and 2sig
			//events that are +/-2 sigma from expected p+8Be(2+) IM2
			//std::cout << "entering p8Be2plus_sig2_check" << std::endl;
			//first, fill the cumulative histograms:
			hDalitzInvMass_p8Be2_2->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_p8Be2_2->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_p8Be2_2->Fill(*ExSPS, *SABREsumE);

			hExSPS_p8Be2_2->Fill(*ExSPS);
			hRecEx_p8Be2_2->Fill(*RecEx);
			hRecEx_VS_ExSPS_p8Be2_2->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_p8Be2_2->Fill(*costhetakT);
			hExEt_T_p8Be2_2->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_p8Be2_2->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_p8Be2_2->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_p8Be2_2->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_p8Be2_2->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_p8Be2_2->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);

			hCosThetaK_Y1_p8Be2_2->Fill(*costhetakY1);
			hExEt_Y1_p8Be2_2->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_p8Be2_2->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_p8Be2_2->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_p8Be2_2->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_p8Be2_2->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_p8Be2_2->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_p8Be2_2->Fill(*costhetakY2);
			hExEt_Y2_p8Be2_2->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_p8Be2_2->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_p8Be2_2->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_p8Be2_2->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_p8Be2_2->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_p8Be2_2->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}
		if(p8Be2plus_sig3_check){ //else if would act as exclusive gates, but we want 2sig to include 1sig and 3sig to include 1sig and 2sig
			//events that are +/-3 sigma from expected p+8Be(2+) IM2
			//std::cout << "entering p8Be2plus_sig3_check" << std::endl;
			//first, fill the cumulative histograms:
			hDalitzInvMass_p8Be2_3->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_p8Be2_3->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_p8Be2_3->Fill(*ExSPS, *SABREsumE);

			hExSPS_p8Be2_3->Fill(*ExSPS);
			hRecEx_p8Be2_3->Fill(*RecEx);
			hRecEx_VS_ExSPS_p8Be2_3->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_p8Be2_3->Fill(*costhetakT);
			hExEt_T_p8Be2_3->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_p8Be2_3->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_p8Be2_3->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_p8Be2_3->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_p8Be2_3->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_p8Be2_3->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);

			hCosThetaK_Y1_p8Be2_3->Fill(*costhetakY1);
			hExEt_Y1_p8Be2_3->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_p8Be2_3->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_p8Be2_3->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_p8Be2_3->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_p8Be2_3->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_p8Be2_3->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_p8Be2_3->Fill(*costhetakY2);
			hExEt_Y2_p8Be2_3->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_p8Be2_3->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_p8Be2_3->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_p8Be2_3->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_p8Be2_3->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_p8Be2_3->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}

		//finally, check a+5Li(gs):
		if(a5Ligs_sig1_check){
			//events that are +/- 1 sigma from expected a+5Li(gs) IM2
			//std::cout << "entering a5Ligs_sig1_check" << std::endl;
			//first, fill the cumulative histograms:
			hDalitzInvMass_a5Ligs_1->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_a5Ligs_1->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_a5Ligs_1->Fill(*ExSPS, *SABREsumE);

			hExSPS_a5Ligs_1->Fill(*ExSPS);
			hRecEx_a5Ligs_1->Fill(*RecEx);
			hRecEx_VS_ExSPS_a5Ligs_1->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_a5Ligs_1->Fill(*costhetakT);
			hExEt_T_a5Ligs_1->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_a5Ligs_1->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_a5Ligs_1->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_a5Ligs_1->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_a5Ligs_1->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_a5Ligs_1->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);

			hCosThetaK_Y1_a5Ligs_1->Fill(*costhetakY1);
			hExEt_Y1_a5Ligs_1->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_a5Ligs_1->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_a5Ligs_1->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_a5Ligs_1->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_a5Ligs_1->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_a5Ligs_1->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_a5Ligs_1->Fill(*costhetakY2);
			hExEt_Y2_a5Ligs_1->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_a5Ligs_1->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_a5Ligs_1->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_a5Ligs_1->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_a5Ligs_1->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_a5Ligs_1->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}
		if(a5Ligs_sig2_check){ //else if would act as exclusive gates, but we want 2sig to include 1sig and 3sig to include 1sig and 2sig
			//events that are +/-2 sigma from expected a+5Li(gs) IM2
			//std::cout << "entering a5Ligs_sig2_check" << std::endl;
			//first, fill the cumulative histograms:
			hDalitzInvMass_a5Ligs_2->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_a5Ligs_2->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_a5Ligs_2->Fill(*ExSPS, *SABREsumE);

			hExSPS_a5Ligs_2->Fill(*ExSPS);
			hRecEx_a5Ligs_2->Fill(*RecEx);
			hRecEx_VS_ExSPS_a5Ligs_2->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_a5Ligs_2->Fill(*costhetakT);
			hExEt_T_a5Ligs_2->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_a5Ligs_2->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_a5Ligs_2->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_a5Ligs_2->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_a5Ligs_2->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_a5Ligs_2->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);

			hCosThetaK_Y1_a5Ligs_2->Fill(*costhetakY1);
			hExEt_Y1_a5Ligs_2->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_a5Ligs_2->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_a5Ligs_2->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_a5Ligs_2->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_a5Ligs_2->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_a5Ligs_2->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_a5Ligs_2->Fill(*costhetakY2);
			hExEt_Y2_a5Ligs_2->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_a5Ligs_2->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_a5Ligs_2->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_a5Ligs_2->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_a5Ligs_2->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_a5Ligs_2->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}
		if(a5Ligs_sig3_check){ //else if would act as exclusive gates, but we want 2sig to include 1sig and 3sig to include 1sig and 2sig
			//events that are +/-3 sigma from expected a+5Li(gs) IM2
			//std::cout << "entering a5Ligs_sig3_check" << std::endl;
			//first, fill the cumulative histograms:
			hDalitzInvMass_a5Ligs_3->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_a5Ligs_3->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_a5Ligs_3->Fill(*ExSPS, *SABREsumE);

			hExSPS_a5Ligs_3->Fill(*ExSPS);
			hRecEx_a5Ligs_3->Fill(*RecEx);
			hRecEx_VS_ExSPS_a5Ligs_3->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_a5Ligs_3->Fill(*costhetakT);
			hExEt_T_a5Ligs_3->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_a5Ligs_3->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_a5Ligs_3->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_a5Ligs_3->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_a5Ligs_3->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_a5Ligs_3->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);

			hCosThetaK_Y1_a5Ligs_3->Fill(*costhetakY1);
			hExEt_Y1_a5Ligs_3->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_a5Ligs_3->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_a5Ligs_3->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_a5Ligs_3->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_a5Ligs_3->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_a5Ligs_3->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_a5Ligs_3->Fill(*costhetakY2);
			hExEt_Y2_a5Ligs_3->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_a5Ligs_3->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_a5Ligs_3->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_a5Ligs_3->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_a5Ligs_3->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_a5Ligs_3->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}
			//inconcistent events:
		if(!p8Be0plus_sig1_check && !p8Be0plus_sig2_check && !p8Be0plus_sig3_check &&
		   !p8Be2plus_sig1_check && !p8Be2plus_sig2_check && !p8Be2plus_sig3_check &&
		   !a5Ligs_sig1_check    && !a5Ligs_sig2_check    && !a5Ligs_sig3_check){
		   	//events that are not +/-3 sigma from any expected channel IM2
			//std::cout << "entering a5Ligs_sig3_check" << std::endl;
			//first, fill the cumulative histograms:
			hDalitzInvMass_none->Fill(*m2_aa, *m2_pa1);
			hDalitzInvMass_none->Fill(*m2_aa, *m2_pa2);

			hSABREsumE_vs_ExSPS_none->Fill(*ExSPS, *SABREsumE);

			hExSPS_none->Fill(*ExSPS);
			hRecEx_none->Fill(*RecEx);
			hRecEx_VS_ExSPS_none->Fill(*ExSPS, *RecEx);
		
			hCosThetaK_T_none->Fill(*costhetakT);
			hExEt_T_none->Fill((*E_xT)/(*E_xT + *E_yT));
			hEyEt_T_none->Fill((*E_yT)/(*E_xT + *E_yT));
			hExEt_CosThetaK_T_none->Fill(*costhetakT, (*E_xT)/(*E_xT + *E_yT));
			hEyEt_CosThetaK_T_none->Fill(*costhetakT, (*E_yT)/(*E_xT + *E_yT));
			hExEt_Et_T_none->Fill((*E_xT + *E_yT), (*E_xT)/(*E_xT + *E_yT));
			hExSPS_ExEt_T_none->Fill((*E_xT)/(*E_xT + *E_yT), *ExSPS);

			hCosThetaK_Y1_none->Fill(*costhetakY1);
			hExEt_Y1_none->Fill((*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_Y1_none->Fill((*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_CosThetaK_Y1_none->Fill(*costhetakY1, (*E_xY1)/(*E_xY1 + *E_yY1));
			hEyEt_CosThetaK_Y1_none->Fill(*costhetakY1, (*E_yY1)/(*E_xY1 + *E_yY1));
			hExEt_Et_Y1_none->Fill((*E_xY1 + *E_yY1), (*E_xY1)/(*E_xY1 + *E_yY1));
			hExSPS_ExEt_Y1_none->Fill((*E_xY1)/(*E_xY1 + *E_yY1), *ExSPS);

			hCosThetaK_Y2_none->Fill(*costhetakY2);
			hExEt_Y2_none->Fill((*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_Y2_none->Fill((*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_CosThetaK_Y2_none->Fill(*costhetakY2, (*E_xY2)/(*E_xY2 + *E_yY2));
			hEyEt_CosThetaK_Y2_none->Fill(*costhetakY2, (*E_yY2)/(*E_xY2 + *E_yY2));
			hExEt_Et_Y2_none->Fill((*E_xY2 + *E_yY2), (*E_xY2)/(*E_xY2 + *E_yY2));
			hExSPS_ExEt_Y2_none->Fill((*E_xY2)/(*E_xY2 + *E_yY2), *ExSPS);
		}

	}
	//std::cout << "out of loop" << std::endl;

	outfile->cd();
	outfile->Write();
	outfile->Close();
	delete outfile;
	infile->Close();
	delete infile;
	std::cout << "Output saved to " << outfilename << std::endl;
}
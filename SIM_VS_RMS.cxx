/*
	Data file name example:			/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_0plus_output.root
	  Histograms contained:			1par/0plus/hSABRE_SABRE3_Ring8ESummary_0plus		(TH1D)

	 Sim file name example:			/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He_6Li3562keV_7500keV.root
	  Histograms contained:			SABRE/SABRE3/Summary/hSABRE3_Ring8Summary



	Need to ensure same binning in both histograms:
			sim histo has 0-20 MeV w/ 2000 bins --> 10 keV/bin
			data histo has 0-20 MeV w/ 2000 bins --> 10 keV/bin

	This resolution may need to change and we can always force simulation to match what is best for data.


	Write individual functions here to handle specific histogram cases.
*/

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TString.h>

void Lithium6_1plus(){
	//uncomment below for surface laptop
	// TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = "1par/1plus/hSABRE_SABRE3_Ring8ESummary_1plus";

	// TString simFilePath = "/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He_6Li3562keV_7500keV.root";
	// TString simHistLocalPath = "SABRE/SABRE3/Summary/hSABRE3_Ring8Summary";

	//uncomment below for DESKTOP
	TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = "1par/1plus/hSABRE_SABRE3_Ring7ESummary_1plus";

	TString simFilePath = "/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He_6Ligs_7500keV_ELoss2.root";
	TString simHistLocalPath = "SABRE/SABRE3/Summary/hSABRE3_Ring7Summary";

	//open data file and retrieve histo

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH1 *hData = dynamic_cast<TH1*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram" << std::endl;
		datafile->Close();
		return;
	}
	hData->SetDirectory(0);
	datafile->Close();

	//open sim file and retrieve histo

	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return;
	}

	TH1 *hSim = dynamic_cast<TH1*>(simfile->Get(simHistLocalPath));
	if(!hSim){
		std::cerr << "Error retrieving sim histogram" << std::endl;
		simfile->Close();
		return;
	}
	hSim->SetDirectory(0);
	simfile->Close();

	//check if binning matches
	if(hData->GetNbinsX() != hSim->GetNbinsX()){
		std::cerr << "Histogram bin counts do not match" << std::endl;
		return;
	}


	//prep output file for scaling factors
	// std::ofstream outfile("Lithium6_1plus_scales.txt");
	// outfile << "# Bin\tCenter\tData\tSim\tScaleFactor" << std::endl;


	//scale sim bin-by-bin to match data
	// int nBins = hData->GetNbinsX();
	// for(int i=1; i <= nBins; i++){
	// 	double dataContent = hData->GetBinContent(i);
	// 	double simContent = hSim->GetBinContent(i);
	// 	double binCenter = hData->GetBinCenter(i);
	// 	double simError = hSim->GetBinError(i);

	// 	double scaleFactor = (simContent != 0) ? (dataContent/simContent) : 0.;
	// 	if(simContent != 0){
	// 		hSim->SetBinContent(i, simContent*scaleFactor);
	// 		hSim->SetBinError(i, simError*scaleFactor);
	// 	}
	// 	//outfile << i << "\t" << binCenter << "\t" << dataContent << "\t" << simContent << "\t" << scaleFactor << std::endl;
	// }

	//global scale factor to scale sim down to exp
	double dataIntegral = hData->Integral();
	double simIntegral = hSim->Integral();

	double scaleFactor = 0.;

	if(simIntegral > 0){
		scaleFactor = dataIntegral/simIntegral;
		hSim->Scale(scaleFactor);
		std::cout << "Applied global scale factor to sim: " << scaleFactor << std::endl;
	} else {
		std::cerr << "sim histogram has zero integral and thus cannot scale" << std::endl;
		return;
	}


	//outfile.close();

	//draw and overlay
	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);

	const int rebinfactor = 4;

	double maxdata = hData->GetMaximum();
	double maxsim = hSim->GetMaximum();

	double ymax = std::max(maxdata,maxsim)*1.1;
	hData->SetMaximum(ymax);

	hData->SetLineColor(kViolet);
	hData->SetLineWidth(4);
	hData->SetTitle("^{6}Li (1^{+} GS) Data vs Sim;Energy (MeV);Counts");
	std::cout << "fitting data with gaussian:" << std::endl;
	TH1D* hData_rebin = dynamic_cast<TH1D*>(hData->Rebin(rebinfactor));
	hData_rebin->Draw("HIST");
	hData_rebin->Fit("gaus", "R SAME");


	std::cout << "\n\nfitting sim with gaussian:" << std::endl;
	hSim->SetLineColor(kOrange);
	hSim->SetLineWidth(2);
	TH1D* hSim_rebin = dynamic_cast<TH1D*>(hSim->Rebin(rebinfactor));
	hSim_rebin->Draw("HIST SAME");
	hSim_rebin->Fit("gaus", "R SAME");

	//add legend
	TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hData,"Data","l");
	TString simlabel = Form("Sim (scaled x %.3f)",scaleFactor);
	legend->AddEntry(hSim,simlabel,"l");
	legend->Draw();

	c1->Update();

	//save to new root file:
	TString outfile_root_name = "Lithium6_1plus_simVSdata_ELoss2.root";
	TFile *outfile_root = new TFile(outfile_root_name,"RECREATE");
	if(!outfile_root || outfile_root->IsZombie()){
		std::cerr << "Error creating output file" << std::endl;
		return;
	}

	hData->Write("hData_original");
	hSim->Write("hSim_scaled");
	c1->Write("overlay_canvas");
	outfile_root->Close();

	std::cout << "Finished! Output written to: " << outfile_root_name << std::endl;
	std::cout << "Scale factors saved to Lithium6_1plus_scales.txt" << std::endl;
}

void Lithium6_0plus(){
//uncomment below for surface laptop
	// TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_0plus_output.root";
	// TString dataHistLocalPath = "1par/0plus/hSABRE_SABRE3_Ring8ESummary_0plus";

	// TString simFilePath = "/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He_6Li3562keV_7500keV.root";
	// TString simHistLocalPath = "SABRE/SABRE3/Summary/hSABRE3_Ring8Summary";

	//uncomment below for DESKTOP
	TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_0plus_output.root";
	TString dataHistLocalPath = "1par/0plus/hSABRE_SABRE3_Ring9ESummary_0plus";

	TString simFilePath = "/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He_6Li3562keV_7500keV_ELoss2.root";
	TString simHistLocalPath = "SABRE/SABRE3/Summary/hSABRE3_Ring9Summary";

	//open data file and retrieve histo

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH1 *hData = dynamic_cast<TH1*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram" << std::endl;
		datafile->Close();
		return;
	}
	hData->SetDirectory(0);
	datafile->Close();

	//open sim file and retrieve histo

	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return;
	}

	TH1 *hSim = dynamic_cast<TH1*>(simfile->Get(simHistLocalPath));
	if(!hSim){
		std::cerr << "Error retrieving sim histogram" << std::endl;
		simfile->Close();
		return;
	}
	hSim->SetDirectory(0);
	simfile->Close();

	//check if binning matches
	if(hData->GetNbinsX() != hSim->GetNbinsX()){
		std::cerr << "Histogram bin counts do not match" << std::endl;
		return;
	}


	//prep output file for scaling factors
	// std::ofstream outfile("Lithium6_0plus_scales.txt");
	// outfile << "# Bin\tCenter\tData\tSim\tScaleFactor" << std::endl;


	//scale sim bin-by-bin to match data
	// int nBins = hData->GetNbinsX();
	// for(int i=1; i <= nBins; i++){
	// 	double dataContent = hData->GetBinContent(i);
	// 	double simContent = hSim->GetBinContent(i);
	// 	double binCenter = hData->GetBinCenter(i);
	// 	double simError = hSim->GetBinError(i);

	// 	double scaleFactor = (simContent != 0) ? (dataContent/simContent) : 0.;
	// 	if(simContent != 0){
	// 		hSim->SetBinContent(i, simContent*scaleFactor);
	// 		hSim->SetBinError(i, simError*scaleFactor);
	// 	}
	// 	//outfile << i << "\t" << binCenter << "\t" << dataContent << "\t" << simContent << "\t" << scaleFactor << std::endl;
	// }

	//global scale factor to scale sim down to exp
	double dataIntegral = hData->Integral();
	double simIntegral = hSim->Integral();

	double scaleFactor = 0.;

	if(simIntegral > 0){
		scaleFactor = dataIntegral/simIntegral;
		hSim->Scale(scaleFactor);
		std::cout << "Applied global scale factor to sim: " << scaleFactor << std::endl;
	} else {
		std::cerr << "sim histogram has zero integral and thus cannot scale" << std::endl;
		return;
	}


	//outfile.close();

	//draw and overlay
	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);

	const int rebinfactor = 4;

	double maxdata = hData->GetMaximum();
	double maxsim = hSim->GetMaximum();

	double ymax = std::max(maxdata,maxsim)*4.1;
	hData->SetMaximum(ymax);

	hData->SetLineColor(kViolet);
	hData->SetLineWidth(4);
	hData->SetTitle("^{6}Li (0^{+} E = 3.562 MeV) Data vs Sim;Energy (MeV);Counts");
	std::cout << "fitting data with gaussian:" << std::endl;
	TH1D* hData_rebin = dynamic_cast<TH1D*>(hData->Rebin(rebinfactor));
	hData_rebin->Draw("HIST");
	hData_rebin->Fit("gaus", "R SAME");


	std::cout << "\n\nfitting sim with gaussian:" << std::endl;
	hSim->SetLineColor(kOrange);
	hSim->SetLineWidth(2);
	TH1D* hSim_rebin = dynamic_cast<TH1D*>(hSim->Rebin(rebinfactor));
	hSim_rebin->Draw("HIST SAME");
	hSim_rebin->Fit("gaus", "R SAME");

	//add legend
	TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hData,"Data","l");
	TString simlabel = Form("Sim (scaled x %.3f)",scaleFactor);
	legend->AddEntry(hSim,simlabel,"l");
	legend->Draw();

	c1->Update();

	//save to new root file:
	TString outfile_root_name = "Lithium6_0plus_simVSdata_ELoss2.root";
	TFile *outfile_root = new TFile(outfile_root_name,"RECREATE");
	if(!outfile_root || outfile_root->IsZombie()){
		std::cerr << "Error creating output file" << std::endl;
		return;
	}

	hData->Write("hData_original");
	hSim->Write("hSim_scaled");
	c1->Write("overlay_canvas");
	outfile_root->Close();

	std::cout << "Finished! Output written to: " << outfile_root_name << std::endl;
	std::cout << "Scale factors saved to Lithium6_0plus_scales.txt" << std::endl;
}
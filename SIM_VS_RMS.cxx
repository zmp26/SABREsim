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
#include <TPaveText.h>
#include <cmath>

void Lithium6_1plus(int ring, TString beamstring){

	TString anglestring = "17282228";

	//uncomment below for surface laptop
	// TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = Form("1par/1plus/hSABRE_SABRE3_Ring%dESummary_1plus",ring);

	// TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%s_histos.root", anglestring.Data(), beamstring.Data());
	// TString simHistLocalPath = Form("SABRE/SABRE3/Summary/hSABRE3_Ring%dSummary",ring);

	//std::cout << simFilePath << std::endl;

	//uncomment below for DESKTOP
	TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = Form("1par/1plus/hSABRE_SABRE3_Ring%dESummary_1plus",ring);

	TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%s_histos.root", anglestring.Data(), beamstring.Data());
	TString simHistLocalPath = Form("SABRE/SABRE3/Summary/hSABRE3_Ring%dSummary",ring);

	//for our sps acceptances, rings 7 and 8 on SABRE3 illuminate
	// if(ring != 7 && ring != 8) {
	// 	std::cout << "ring not found, try 7 or 8" << std::endl;
	// 	return;
	// }

	//prep output file name using ring:
	TString outfile_root_name;
	//outfile_root_name = Form("Lithium6_1plus_simVSdata_ring%d.root",ring);//"Lithium6_0plus_simVSdata_ELoss2.root";
	outfile_root_name = Form("Lithium6_1plus_simVSdata_%s_theta%s_ring%d.root", beamstring.Data(), anglestring.Data(), ring);



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

	double ymax = std::max(maxdata,maxsim)*1.1*rebinfactor;
	std::cout << "maxdata = " << maxdata << "\tmaxsim = " << maxsim << std::endl;
	std::cout << "ymax = " << ymax << std::endl;
	hData->SetMaximum(ymax);

	hData->SetLineColor(kViolet);
	hData->SetLineWidth(4);
	TString title = Form("^{6}Li (1^{+} GS) Data vs Sim (Ring %d);Energy (MeV);Counts",ring);
	hData->SetTitle(title);
	//std::cout << "fitting data with gaussian:" << std::endl;
	TH1D* hData_rebin = dynamic_cast<TH1D*>(hData->Rebin(rebinfactor));
	hData_rebin->GetXaxis()->SetRangeUser(2,4);
	hData_rebin->GetYaxis()->SetRangeUser(0,1100);
	hData_rebin->Draw("HIST");
	//hData_rebin->Fit("gaus", "R SAME");


	//std::cout << "\n\nfitting sim with gaussian:" << std::endl;
	hSim->SetLineColor(kOrange);
	hSim->SetLineWidth(2);
	TH1D* hSim_rebin = dynamic_cast<TH1D*>(hSim->Rebin(rebinfactor));
	hSim_rebin->GetXaxis()->SetRangeUser(2,4);
	hSim_rebin->GetYaxis()->SetRangeUser(0,1100);
	hSim_rebin->Draw("HIST SAME");
	//hSim_rebin->Fit("gaus", "R SAME");

	//calculate chi2:
	double chi2 = hData->Chi2Test(hSim, "CHI2/NDF UW");
	double chi2_rebin = hData_rebin->Chi2Test(hSim_rebin, "CHI2/NDF UW");

	//add legend
	TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hData,"Data","l");
	TString simlabel = Form("Sim (scaled x %.3f)",scaleFactor);
	legend->AddEntry(hSim,simlabel,"l");
	legend->Draw();

	//add text box to differentiate pngs!
	TPaveText *pt = new TPaveText(0.65, 0.63, 0.88, 0.74, "NDC");
	pt->AddText(Form("%s SPSTheta%s", beamstring.Data(), anglestring.Data()));
	pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2));
	//pt->AddText(Form("#chi_{rebin}^{2}/NDF = %.3f",chi2_rebin));
	pt->SetFillColorAlpha(kWhite, 0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	c1->Update();

	//save to new root file:
	//TString outfile_root_name = "Lithium6_1plus_simVSdata_ELoss2.root";
	// TFile *outfile_root = new TFile(outfile_root_name,"RECREATE");
	// if(!outfile_root || outfile_root->IsZombie()){
	// 	std::cerr << "Error creating output file" << std::endl;
	// 	return;
	// }

	//hData->Write("hData_original");
	//hSim->Write("hSim_scaled");
	//c1->Write("overlay_canvas");
	c1->SaveAs(Form("Lithium6_1plus_simVSdata_%s_theta%s_ring%d.png", beamstring.Data(), anglestring.Data(), ring));
	//outfile_root->Close();

	//std::cout << "Finished! Output written to: " << outfile_root_name << std::endl;
	//std::cout << "Scale factors saved to Lithium6_1plus_scales.txt" << std::endl;
}

void Lithium6_1plus_auto(){
	std::vector<int> rings = {7,8,9};
	std::vector<TString> beamstrings = {
										"fixedpoint",
										"gaus001",
										"gaus002",
										"gaus003",
										"gaus004",
										"gaus005"
										// "gaus006",
										// "gaus007",
										// "gaus008",
										// "gaus009",
										// "gaus010"
									};

	for(size_t r=0; r<rings.size(); r++){
		for(size_t b=0; b<beamstrings.size(); b++){
			Lithium6_1plus(rings[r],beamstrings[b]);
		}
	}
}

void Lithium6_1plus_sabrehits(TString beamstring){

	TString anglestring = "17282228";

	//uncomment below for DESKTOP
	TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	TString simHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	//uncomment below for LAPTOP
	// TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	// TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	// TString simHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";



	//prep output file name
	//TString outfile_root_name = Form("Lithium6_1plus_simVSdata_sabrehits_%s.root",beamstring.Data());


	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram" << std::endl;
		datafile->Close();
		return;
	}
	hData->SetDirectory(0);
	datafile->Close();



	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return;
	}

	TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
	if(!hSim){
		std::cerr << "Error retrieving sim histogram" << std::endl;
		simfile->Close();
		return;
	}
	hSim->SetDirectory(0);
	simfile->Close();



	if((hData->GetNbinsX() != hSim->GetNbinsX()) || (hData->GetNbinsY() != hSim->GetNbinsY())){
		std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
		return;
	}



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

	//calculate chi2:
	double chi2 = hData->Chi2Test(hSim, "CHI2/NDF UW");


	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);



	hSim->Draw();

	TPaveText *pt = new TPaveText(0.65, 0.75, 0.88, 0.88, "NDC");
	pt->AddText(Form("%s (scale factor = %.3f)", beamstring.Data(), scaleFactor));
	pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	c1->Update();

	// TFile* outfile_root = new TFile(outfile_root_name,"RECREATE");
	// if(!outfile_root || outfile_root->IsZombie()){
	// 	std::cerr << "Error creating output file" << std::endl;
	// 	return;
	// }

	//hData->Write("hData_original");
	//hSim->Write("hSim_scaled");
	c1->SaveAs(Form("Lithium6_1plus_simVSdata_sabrehits_%s.png",beamstring.Data()));
	c1->Clear();
	hData->Draw();
	pt = new TPaveText(0.65, 0.75, 0.88, 0.88, "NDC");
	pt->AddText("LiFha exp 1par ^{6}Li 1^{+} DATA");
	c1->Update();
	c1->SaveAs(Form("Lithium6_1plus_simVSdata_sabrehits_%s.png","1PAREXPDATA"));
	//outfile_root->Close();

	//std::cout << "Finished! Output written to: " << outfile_root_name << std::endl;

}

void Lithium6_1plus_sabrehits_auto(){

	std::vector<TString> beamstrings = {
										"fixedpoint",
										"gaus001",
										"gaus002",
										"gaus003",
										"gaus004",
										"gaus005"
										// "gaus006",
										// "gaus007",
										// "gaus008",
										// "gaus009",
										// "gaus010"
									};

	for(size_t b=0; b<beamstrings.size(); b++){
		Lithium6_1plus_sabrehits(beamstrings[b]);
	}

}

void Lithium6_1plus_sabre3summaries(TString beamstring){

	TString anglestring = "17282228";

	//uncomment below for DESKTOP
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath_rings = "SABRE/SABRE3/Summary/hSABRE3_ERingSummary";
	// TString dataHistLocalPath_wedges = "SABRE/SABRE3/Summary/hSABRE3_EWedgeSummary";

	// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	// TString simHistLocalPath_rings = "SABRE/SABRE3/Summary/hSABRE3_ERingSummary";
	// TString simHistLocalPath_wedges = "SABRE/SABRE3/Summary/hSABRE3_EWedgeSummary";

	//uncomment below for LAPTOP
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath_rings = "SABRE/SABRE3/Summary/hSABRE3_ERingSummary";
	TString dataHistLocalPath_wedges = "SABRE/SABRE3/Summary/hSABRE3_EWedgeSummary";

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	TString simHistLocalPath_rings = "SABRE/SABRE3/Summary/hSABRE3_ERingSummary";
	TString simHistLocalPath_wedges = "SABRE/SABRE3/Summary/hSABRE3_EWedgeSummary";

	//prep output file name:
	//TString outfile_root_name = Form("Lithium6_1plus_simVSdata_sabre3summaries_%s.root",beamstring.Data());


	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH1* hDataRings = dynamic_cast<TH1*>(datafile->Get(dataHistLocalPath_rings));
	TH1* hDataWedges = dynamic_cast<TH1*>(datafile->Get(dataHistLocalPath_wedges));
	if(!hDataRings || !hDataWedges){
		std::cerr << "Error retrieving data histograms!" << std::endl;
		datafile->Close();
		return;
	}
	hDataRings->SetDirectory(0);
	hDataWedges->SetDirectory(0);
	datafile->Close();


	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return;
	}

	TH1 *hSimRings = dynamic_cast<TH1*>(simfile->Get(simHistLocalPath_rings));
	TH1 *hSimWedges = dynamic_cast<TH1*>(simfile->Get(simHistLocalPath_wedges));
	if(!hSimRings || !hSimWedges){
		std::cerr << "Error retrieving sim histograms!" << std::endl;
		simfile->Close();
		return;
	}
	hSimRings->SetDirectory(0);
	hSimWedges->SetDirectory(0);
	simfile->Close();


	if( (hDataRings->GetNbinsX() != hSimRings->GetNbinsX()) || (hDataWedges->GetNbinsX() != hSimWedges->GetNbinsX()) ){
		std::cerr << "Histogram binning does not match on rings and/or wedges!" << std::endl;
		TString temp = Form("hDataRings = %d\thSimRings = %d\nhDataWedges = %d\thSimWedges = %d\n\n", hDataRings->GetNbinsX(), hSimRings->GetNbinsX(), hDataWedges->GetNbinsX(), hSimWedges->GetNbinsX());
		std::cerr << temp.Data();
		return;
	}


	double dataIntegralRings = hDataRings->Integral();
	double dataIntegralWedges = hDataWedges->Integral();
	double simIntegralRings = hSimRings->Integral();
	double simIntegralWedges = hSimWedges->Integral();

	double scaleFactorRings = 0.;
	double scaleFactorWedges = 0.;

	if(simIntegralRings > 0 && simIntegralWedges > 0){
		scaleFactorRings = dataIntegralRings/simIntegralRings;
		scaleFactorWedges = dataIntegralWedges/simIntegralWedges;
		hSimRings->Scale(scaleFactorRings);
		hSimWedges->Scale(scaleFactorWedges);
		std::cout << "applied global scale factors\trings = " << scaleFactorRings << "\twedges = " << scaleFactorWedges << std::endl;
	} else {
		std::cerr << "sim histogram for rings and/or wedges has zero integral and thus cannot be scaled" << std::endl;
		return;
	}

	//prep output root file
	// TFile *outfile_root = new TFile(outfile_root_name, "RECREATE");
	// if(!outfile_root || outfile_root->IsZombie()){
	// 	std::cerr << "Error creating output file" << std::endl;
	// 	return;
	// }

	const int rebinfactor = 4;

	double maxdatarings = hDataRings->GetMaximum();
	double maxdatawedges = hDataWedges->GetMaximum();

	double maxsimrings = hSimRings->GetMaximum();
	double maxsimwedges = hSimWedges->GetMaximum();

	double ymaxrings = std::max(maxdatarings,maxsimrings)*1.1*rebinfactor;
	double ymaxwedges = std::max(maxdatawedges,maxsimwedges)*1.1*rebinfactor;

	hDataRings->SetMaximum(ymaxrings);
	hDataWedges->SetMaximum(ymaxwedges);

	//start drawing w/ canvas

	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);

	//calculate chi squareds here:
	double chi2_rings = hDataRings->Chi2Test(hSimRings, "CHI2/NDF UW");
	double chi2_wedges = hDataWedges->Chi2Test(hSimWedges, "CHI2/NDF UW");

	hDataRings->SetLineColor(kViolet);
	hDataRings->SetLineWidth(4);
	TString title = "^{6}Li (1^{+} GS) Data vs Sim (SABRE3 Ring Summary);Energy (MeV)";
	hDataRings->SetTitle(title);
	TH1D* hDataRings_rebin = dynamic_cast<TH1D*>(hDataRings->Rebin(rebinfactor));
	hDataRings_rebin->GetXaxis()->SetRangeUser(2,4);
	hDataRings_rebin->GetYaxis()->SetRangeUser(0,1400);
	hDataRings_rebin->Draw("HIST");

	hSimRings->SetLineColor(kOrange);
	hSimRings->SetLineWidth(2);
	TH1D* hSimRings_rebin = dynamic_cast<TH1D*>(hSimRings->Rebin(rebinfactor));
	hSimRings_rebin->GetXaxis()->SetRangeUser(2,4);
	hSimRings_rebin->GetYaxis()->SetRangeUser(0,1400);
	hSimRings_rebin->Draw("HIST SAME");


	//add legend
	TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hDataRings,"Data","l");
	TString simringlabel = Form("Sim (scaled x %.3f)",scaleFactorRings);
	legend->AddEntry(hSimRings,simringlabel,"l");
	legend->Draw();

	//add text box to differentiate pngs!
	TPaveText *pt = new TPaveText(0.65,0.63,0.88,0.74,"NDC");
	pt->AddText(beamstring);
	pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2_rings));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	c1->Update();
	//save to png file:
	c1->SaveAs(Form("Lithium6_1plus_simVSdata_SABRE3_RingSummary_%s.png",beamstring.Data()));



	//clear canvas, do same thing but for wedges this time!
	c1->Clear();

	hDataWedges->SetLineColor(kViolet);
	hDataWedges->SetLineWidth(4);
	title = "^{6}Li (1^{+} GS) Data vs Sim (SABRE3 Wedge Summary);Energy (MeV)";
	hDataWedges->SetTitle(title);
	TH1D* hDataWedges_rebin = dynamic_cast<TH1D*>(hDataWedges->Rebin(rebinfactor));
	hDataWedges_rebin->GetXaxis()->SetRangeUser(2,4);
	hDataWedges_rebin->GetYaxis()->SetRangeUser(0,1400);
	hDataWedges_rebin->Draw("HIST");

	hSimWedges->SetLineColor(kOrange);
	hSimWedges->SetLineWidth(2);
	TH1D* hSimWedges_rebin = dynamic_cast<TH1D*>(hSimWedges->Rebin(rebinfactor));
	hSimWedges_rebin->GetXaxis()->SetRangeUser(2,4);
	hSimWedges_rebin->GetYaxis()->SetRangeUser(0,1400);
	hSimWedges_rebin->Draw("HIST SAME");

	//add legend:
	legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hDataWedges,"Data","l");
	simringlabel = Form("Sim (scaled x %.3f)",scaleFactorWedges);
	legend->AddEntry(hSimWedges,simringlabel,"l");
	legend->Draw();

	//add text box to differentiate pngs!
	pt = new TPaveText(0.65,0.63,0.88,0.74,"NDC");
	pt->AddText(beamstring);
	pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2_wedges));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	c1->Update();
	//save to png file:
	c1->SaveAs(Form("Lithium6_1plus_simVSdata_SABRE3_WedgeSummary_%s.png",beamstring.Data()));
}

void Lithium6_1plus_sabre3summaries_auto(){
	std::vector<TString> beamstrings = {
										"fixedpoint",
										"gaus001",
										"gaus002",
										"gaus003",
										"gaus004",
										"gaus005",
										// "gaus006",
										// "gaus007",
										// "gaus008",
										// "gaus009",
										// "gaus010"
									};

	for(const auto& bs : beamstrings){
		Lithium6_1plus_sabre3summaries(bs);
	}
}

void Lithium6_1plus_pixelhistos(int ringChan, int wedgeChan, TString beamstring){

	int SABRE_ID = 3;

	TString pixelhistoname = Form("hSABRE%d_pixel_r%dw%d_ESummary",SABRE_ID,ringChan,wedgeChan);
	TString anglestring = "17282228";

	//uncomment below for DESKTOP:
	//TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	//TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());

	//uncomment below for LAPTOP:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = Form("SABRE/SABRE%d/Pixels/",SABRE_ID);
	dataHistLocalPath = dataHistLocalPath + pixelhistoname;

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	TString simHistLocalPath = Form("SABRE/SABRE%d/Pixels/",SABRE_ID);
	simHistLocalPath = simHistLocalPath + pixelhistoname;


	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		std::cerr << "dataFilePath = " << dataFilePath << "\n\n";
		return;
	}

	TH1 *hData = dynamic_cast<TH1*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram" << std::endl;
		std::cerr << "dataFilePath = " << dataFilePath << "\n";
		std::cerr << "dataHistLocalPath = " << dataHistLocalPath << "\n\n";
		datafile->Close();
		return;
	}
	hData->SetDirectory(0);
	datafile->Close();



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



	if((hData->GetNbinsX() != hSim->GetNbinsX()) || (hData->GetNbinsY() != hSim->GetNbinsY())){
		std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
		return;
	}


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

	const int rebinfactor = 4;

	double maxdata = hData->GetMaximum();
	double maxsim = hSim->GetMaximum();

	double ymax = std::max(maxdata,maxsim);

	//calculate chi2:
	double chi2 = hData->Chi2Test(hSim, "CHI2/NDF UW");


	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);

	hData->SetLineColor(kViolet);
	hData->SetLineWidth(4);
	hData->SetTitle(pixelhistoname);
	TH1D *hDataRebin = dynamic_cast<TH1D*>(hData->Rebin(rebinfactor));
	hDataRebin->GetXaxis()->SetRangeUser(2,4);
	hDataRebin->GetYaxis()->SetRangeUser(0,1000);
	hDataRebin->Draw("HIST");

	hSim->SetLineColor(kOrange);
	hSim->SetLineWidth(2);
	TH1D* hSimRebin = dynamic_cast<TH1D*>(hSim->Rebin(rebinfactor));
	hSimRebin->GetXaxis()->SetRangeUser(2,4);
	hSimRebin->GetYaxis()->SetRangeUser(0,1000);
	hSimRebin->Draw("HIST SAME");

	TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hDataRebin,"Data","l");
	TString simlabel = Form("Sim (scaled x %.3f)",scaleFactor);
	legend->AddEntry(hSimRebin,simlabel,"l");
	legend->Draw();

	TPaveText *pt = new TPaveText(0.65, 0.63, 0.88, 0.74, "NDC");
	pt->AddText(Form("%s (scale factor = %.3f)", beamstring.Data(), scaleFactor));
	pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	c1->Update();

	c1->SaveAs(Form("Lithium6_1plus_simVSdata_pixelhisto_r%dw%d_theta%s_%s.png",ringChan,wedgeChan,anglestring.Data(),beamstring.Data()));

	delete c1;
	delete hData;
	//delete hDataRebin;
	delete hSim;
	//delete hSimRebin;
}

void Lithium6_1plus_pixelhistos_auto(){
	std::vector<TString> beamstrings = {
										"fixedpoint",
										"gaus001",
										"gaus002",
										"gaus003",
										"gaus004",
										"gaus005",
										// "gaus006",
										// "gaus007",
										// "gaus008",
										// "gaus009",
										// "gaus010"
									};

	std::vector<int> ringChans = {
									71,
									72,
									73	
									};

	std::vector<int> wedgeChans = {
									29,
									30
									};

	for(const auto& bs : beamstrings){
		for(const auto& ring : ringChans){
			for(const auto& wedge : wedgeChans){
				Lithium6_1plus_pixelhistos(ring,wedge,bs);
			}
		}
	}
}

void Lithium6_1plus_fourpixelchisquared(){

	TString anglestring = "17282228";

	std::vector<TString> beamstrings = {
											"fixed",
											"gaus001",
											"gaus002",
											"gaus003",
											"gaus004",
											"gaus005"
										};

	std::vector<TString> filenames;

	for(const auto& bsx : beamstrings){
		
		for(const auto& bsy : beamstrings){
			
			TString filename = Form("kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%sx_%sy_histos.root", anglestring.Data(), bsx.Data(), bsy.Data());
			filenames.push_back(filename);

		}

	}

	//---------------------------------------------------
	//				Establish data values
	//---------------------------------------------------

	//uncomment for DESKTOP:
	TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";

	TString path_pix_r71_w29 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r71_w29";
	TString path_pix_r72_w29 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r72_w29";
	TString path_pix_r72_w30 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r72_w30";
	TString path_pix_r71_w30 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r71_w30";

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH1 *hpix_r71_w29 = dynamic_cast<TH1*>(datafile->Get(path_pix_r71_w29));
	TH1 *hpix_r72_w29 = dynamic_cast<TH1*>(datafile->Get(path_pix_r72_w29));
	TH1 *hpix_r72_w30 = dynamic_cast<TH1*>(datafile->Get(path_pix_r72_w30));
	TH1 *hpix_r71_w30 = dynamic_cast<TH1*>(datafile->Get(path_pix_r71_w30));
	if(!hpix_r71_w29 || !hpix_r72_w29 || !hpix_r72_w30 || !hpix_r71_w30){
		std::cerr << "Error retrieving at least one data histogram!" << std::endl;
		return;
	}

	// hpix_r71_w29->SetDirectory(0);
	// hpix_r72_w29->SetDirectory(0);
	// hpix_r72_w30->SetDirectory(0);
	// hpix_r71_w30->SetDirectory(0);
	datafile->Close();

	double counts_pix_r71_w29 = hpix_r71_w29->GetEntries();
	double counts_pix_r72_w29 = hpix_r72_w29->GetEntries();
	double counts_pix_r72_w30 = hpix_r72_w30->GetEntries();
	double counts_pix_r71_w30 = hpix_r71_w30->GetEntries();

	double fourpixsum = counts_pix_r71_w29 + counts_pix_r72_w29 + counts_pix_r72_w30 + counts_pix_r71_w30;

	std::vector<double> fourpix_relcounts = {counts_pix_r71_w29/fourpixsum, counts_pix_r72_w29/fourpixsum, counts_pix_r72_w30/fourpixsum, counts_pix_r71_w30/fourpixsum};

	for(const auto& fn : filenames){

		//establish sim values:
		TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/%s",fn.Data());

		TFile *simfile = new TFile(simFilePath,"READ");
		if(!simfile || simfile->IsZombie()){
			std::cerr << "Error opening sim file " << simFilePath << "\n\n";
			continue;
		}

		TH1 *hpix_r71_w29_sim = dynamic_cast<TH1*>(simfile->GetPath(path_pix_r71_w29));
		TH1 *hpix_r72_w29_sim = dynamic_cast<TH1*>(simfile->GetPath(path_pix_r72_w29));
		TH1 *hpix_r72_w30_sim = dynamic_cast<TH1*>(simfile->GetPath(path_pix_r72_w30));
		TH1 *hpix_r71_w30_sim = dynamic_cast<TH1*>(simfile->GetPath(path_pix_r71_w30));
		if(!hpix_r71_w29_sim || !hpix_r72_w29_sim || !hpix_r72_w30_sim || !hpix_r71_w30_sim){
			std::cerr << "Error retrieving at least one sim histogram!" << std::endl;
			continue;
		}

		// hpix_r71_w29_sim->SetDirectory(0);
		// hpix_r72_w29_sim->SetDirectory(0);
		// hpix_r72_w30_sim->SetDirectory(0);
		// hpix_r71_w30_sim->SetDirectory(0);
		simfile->Close();

		if( (hpix_r71_w29->GetNbinsX() != hpix_r71_w29_sim->GetNbinsX()) || (hpix_r72_w29->GetNbinsX() != hpix_r72_w29_sim->GetNbinsX()) || (hpix_r72_w30->GetNbinsX() != hpix_r72_w30_sim->GetNbinsX()) || (hpix_r71_w30->GetNbinsX() != hpix_r71_w30_sim->GetNbinsX())){
			std::cerr << "Histogram binning does not match for at least one data/sim histogram pair\n";
			continue;
		}

		double simIntegral_r71_w29 = hpix_r71_w29_sim->Integral();
		double simIntegral_r72_w29 = hpix_r72_w29_sim->Integral();
		double simIntegral_r72_w30 = hpix_r72_w30_sim->Integral();
		double simIntegral_r71_w30 = hpix_r71_w30_sim->Integral();

		double scaleFactor_r71_w29;// = counts_pix_r71_w29/simIntegral_r71_w29;
		double scaleFactor_r72_w29;// = counts_pix_r72_w29/simIntegral_r72_w29;
		double scaleFactor_r72_w30;// = counts_pix_r72_w30/simIntegral_r72_w30;
		double scaleFactor_r71_w30;// = counts_pix_r71_w30/simIntegral_r71_w30;

		if(simIntegral_r71_w29!=0 && simIntegral_r72_w29!=0 && simIntegral_r72_w30!=0 && simIntegral_r71_w30!=0){

			scaleFactor_r71_w29 = counts_pix_r71_w29/simIntegral_r71_w29;
			scaleFactor_r72_w29 = counts_pix_r72_w29/simIntegral_r72_w29;
			scaleFactor_r72_w30 = counts_pix_r72_w30/simIntegral_r72_w30;
			scaleFactor_r71_w30 = counts_pix_r71_w30/simIntegral_r71_w30;

			hpix_r71_w29_sim->Scale(scaleFactor_r71_w29);
			hpix_r72_w29_sim->Scale(scaleFactor_r72_w29);
			hpix_r72_w30_sim->Scale(scaleFactor_r72_w30);
			hpix_r71_w30_sim->Scale(scaleFactor_r71_w30);

		} else {
			std::cerr << "sim histogram has zero integral, cannot be scaled...continuing...\n";
			continue;
		}


		double counts_sim_pix_r71_w29 = hpix_r71_w29_sim->Integral();
		double counts_sim_pix_r72_w29 = hpix_r72_w29_sim->Integral();
		double counts_sim_pix_r72_w30 = hpix_r72_w30_sim->Integral();
		double counts_sim_pix_r71_w30 = hpix_r71_w30_sim->Integral();

		double fourpixsum_sim = counts_sim_pix_r71_w29 + counts_sim_pix_r72_w29 + counts_sim_pix_r72_w30 + counts_sim_pix_r71_w30;

		std::vector<double> fourpix_relcounts_sim = {counts_sim_pix_r71_w29/fourpixsum_sim, counts_sim_pix_r72_w29/fourpixsum_sim, counts_sim_pix_r72_w30/fourpixsum_sim, counts_sim_pix_r71_w30/fourpixsum_sim};

		//calculate the chi squared for this sim file:

		double chi2 = (std::pow((counts_pix_r71_w29-counts_sim_pix_r71_w29),2)/std::pow((),2)) + (std::pow((),2)/) + (std::pow((),2)/) + (std::pow((),2)/);

	}

}

void Lithium6_0plus(int ring, TString beamstring){

	TString anglestring = "17282228";

	//uncomment below for surface laptop
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_0plus_output.root";
	TString dataHistLocalPath = Form("1par/0plus/hSABRE_SABRE3_Ring%dESummary_0plus",ring);

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Li3563_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	TString simHistLocalPath = Form("SABRE/SABRE3/Summary/hSABRE3_Ring%dSummary",ring);

	// if(ring != 8 && ring != 9){
	// 	std::cout << "ring histo not there, try 8 or 9" << std::endl;
	// 	return;
	// }

	//prep output file name using ring:
	TString outfile_root_name = Form("Lithium6_0plus_sim%sVSdata_ring%d.root",anglestring.Data(),ring);//"Lithium6_0plus_simVSdata_ELoss2.root";

	//uncomment below for DESKTOP
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_0plus_output.root";
	// TString dataHistLocalPath = Form("1par/0plus/hSABRE_SABRE3_Ring%dESummary_0plus",ring);

	// TString simFilePath = "/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Li3562_7500keV_theta1721_histos.root";
	// TString simHistLocalPath = Form("SABRE/SABRE3/Summary/hSABRE3_Ring%dSummary",ring);

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

	double ymax = std::max(maxdata,maxsim)*rebinfactor*1.1;
	hData->SetMaximum(ymax);

	hData->SetLineColor(kViolet);
	hData->SetLineWidth(4);
	hData->SetTitle(Form("^{6}Li (0^{+} E = 3.563 MeV) Data vs Sim (Ring %d);Energy (MeV);Counts",ring));
	//std::cout << "fitting data with gaussian:" << std::endl;
	TH1D* hData_rebin = dynamic_cast<TH1D*>(hData->Rebin(rebinfactor));
	hData_rebin->GetXaxis()->SetRangeUser(0,4);
	hData_rebin->GetYaxis()->SetRangeUser(0,250);
	hData_rebin->Draw("HIST");
	//hData_rebin->Fit("gaus", "R SAME");


	//std::cout << "\n\nfitting sim with gaussian:" << std::endl;
	hSim->SetLineColor(kOrange);
	hSim->SetLineWidth(2);
	TH1D* hSim_rebin = dynamic_cast<TH1D*>(hSim->Rebin(rebinfactor));
	hSim_rebin->GetXaxis()->SetRangeUser(0,4);
	hSim_rebin->GetYaxis()->SetRangeUser(0,250);
	hSim_rebin->Draw("HIST SAME");
	//hSim_rebin->Fit("gaus", "R SAME");


	//calculate chi2:
	double chi2 = hData->Chi2Test(hSim, "CHI2/NDF UW");
	double chi2_rebin = hData_rebin->Chi2Test(hSim_rebin, "CHI2/NDF UW");

	//add legend
	TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hData,"Data","l");
	TString simlabel = Form("Sim (scaled x %.3f)",scaleFactor);
	legend->AddEntry(hSim,simlabel,"l");
	legend->Draw();

	//add txt box:
	TPaveText *pt = new TPaveText(0.65, 0.63, 0.88, 0.74, "NDC");
	pt->AddText(Form("%s SPSTheta%s", beamstring.Data(), anglestring.Data()));
	pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2));
	//pt->AddText(Form("#chi_{rebin}^{2}/NDF = %.3f",chi2_rebin));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	c1->Update();

	c1->SaveAs(Form("Lithium6_0plus_simVSdata_%s_theta%s_ring%d.png",beamstring.Data(),anglestring.Data(),ring));

	//save to new root file:
	//TString outfile_root_name = "Lithium6_0plus_simVSdata_ELoss2.root";
	// TFile *outfile_root = new TFile(outfile_root_name,"RECREATE");
	// if(!outfile_root || outfile_root->IsZombie()){
	// 	std::cerr << "Error creating output file" << std::endl;
	// 	return;
	// }

	// hData->Write("hData_original");
	// hSim->Write("hSim_scaled");
	// c1->Write("overlay_canvas");
	//outfile_root->Close();

	std::cout << "Finished " << beamstring << " " << ring << "!" << std::endl;
	//std::cout << "Scale factors saved to Lithium6_0plus_scales.txt" << std::endl;
}

void Lithium6_0plus_auto(){
	std::vector<int> rings = {7,8,9};
	std::vector<TString> beamstrings = {
										"fixedpoint",
										"gaus001",
										"gaus002",
										"gaus003",
										"gaus004",
										"gaus005",
										"gaus006",
										"gaus007",
										"gaus008",
										"gaus009",
										"gaus010"
									};

	for(const auto& ring : rings){
		for(const auto& bs : beamstrings){
			Lithium6_0plus(ring,bs);
		}
	}
}

void Lithium6_0plus_sabrehits(TString beamstring){
	//uncomment below for DESKTOP
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_0plus_output.root";
	// TString dataHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Li3563_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	// TString simHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	TString anglestring = "17282228";

	//uncomment below for LAPTOP
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_0plus_output.root";
	TString dataHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Li3563_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	TString simHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";


	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram" << std::endl;
		datafile->Close();
		return;
	}
	hData->SetDirectory(0);
	datafile->Close();


	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return;
	}

	TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
	if(!hSim){
		std::cerr << "Error retrieving sim histogram" << std::endl;
		simfile->Close();
		return;
	}
	hSim->SetDirectory(0);
	simfile->Close();


	if((hData->GetNbinsX() != hSim->GetNbinsX()) || (hData->GetNbinsY() != hSim->GetNbinsY())){
		std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
		return;
	}



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

	//calculate chi2:
	double chi2 = hData->Chi2Test(hSim, "CHI2/NDF UW");


	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);




	hSim->Draw();

	TPaveText *pt = new TPaveText(0.65, 0.75, 0.88, 0.88, "NDC");
	pt->AddText(Form("%s (scale factor = %.3f)", beamstring.Data(), scaleFactor));
	pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	c1->Update();

	c1->SaveAs(Form("Lithium6_0plus_simVSdata_sabrehits_%s.png",beamstring.Data()));
	c1->Clear();
	hData->Draw();
	pt = new TPaveText(0.65, 0.75, 0.88, 0.88, "NDC");
	pt->AddText("LiFha exp 1par ^{6}Li 0^{+} DATA");
	c1->Update();
	c1->SaveAs(Form("Lithium6_0plus_simVSdata_sabrehits_%s.png","1PAREXPDATA"));

}

void Lithium6_0plus_sabrehits_auto(){
	std::vector<TString> beamstrings = {
										"fixedpoint",
										"gaus001",
										"gaus002",
										"gaus003",
										"gaus004",
										"gaus005",
										"gaus006",
										"gaus007",
										"gaus008",
										"gaus009",
										"gaus010"
									};

	for(size_t b=0; b<beamstrings.size(); b++){
		Lithium6_0plus_sabrehits(beamstrings[b]);
	}	
}

void Lithium6_0plus_sabre3summaries(TString beamstring){
	TString anglestring = "17282228";

	//uncomment below for DESKTOP
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_0plus_output.root";
	// TString dataHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Li3563_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	// TString simHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	//uncomment below for LAPTOP
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_0plus_output.root";
	TString dataHistLocalPath_rings = "SABRE/SABRE3/Summary/hSABRE3_ERingSummary";
	TString dataHistLocalPath_wedges = "SABRE/SABRE3/Summary/hSABRE3_EWedgeSummary";

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Li3563_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());
	TString simHistLocalPath_rings = "SABRE/SABRE3/Summary/hSABRE3_ERingSummary";
	TString simHistLocalPath_wedges = "SABRE/SABRE3/Summary/hSABRE3_EWedgeSummary";


	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH1* hDataRings = dynamic_cast<TH1*>(datafile->Get(dataHistLocalPath_rings));
	TH1* hDataWedges = dynamic_cast<TH1*>(datafile->Get(dataHistLocalPath_wedges));
	if(!hDataRings || !hDataWedges){
		std::cerr << "Error retrieving data histograms!" << std::endl;
		datafile->Close();
		return;
	}
	hDataRings->SetDirectory(0);
	hDataWedges->SetDirectory(0);
	datafile->Close();


	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return;
	}

	TH1 *hSimRings = dynamic_cast<TH1*>(simfile->Get(simHistLocalPath_rings));
	TH1 *hSimWedges = dynamic_cast<TH1*>(simfile->Get(simHistLocalPath_wedges));
	if(!hSimRings || !hSimWedges){
		std::cerr << "Error retrieving sim histograms!" << std::endl;
		simfile->Close();
		return;
	}
	hSimRings->SetDirectory(0);
	hSimWedges->SetDirectory(0);
	simfile->Close();


	if( (hDataRings->GetNbinsX() != hSimRings->GetNbinsX()) || (hDataWedges->GetNbinsX() != hSimWedges->GetNbinsX()) ){
		std::cerr << "Histogram binning does not match on rings and/or wedges!" << std::endl;
		TString temp = Form("hDataRings = %d\thSimRings = %d\nhDataWedges = %d\thSimWedges = %d\n\n", hDataRings->GetNbinsX(), hSimRings->GetNbinsX(), hDataWedges->GetNbinsX(), hSimWedges->GetNbinsX());
		std::cerr << temp.Data();
		return;
	}


	double dataIntegralRings = hDataRings->Integral();
	double dataIntegralWedges = hDataWedges->Integral();
	double simIntegralRings = hSimRings->Integral();
	double simIntegralWedges = hSimWedges->Integral();

	double scaleFactorRings = 0.;
	double scaleFactorWedges = 0.;

	if(simIntegralRings > 0 && simIntegralWedges > 0){
		scaleFactorRings = dataIntegralRings/simIntegralRings;
		scaleFactorWedges = dataIntegralWedges/simIntegralWedges;
		hSimRings->Scale(scaleFactorRings);
		hSimWedges->Scale(scaleFactorWedges);
		std::cout << "applied global scale factors\trings = " << scaleFactorRings << "\twedges = " << scaleFactorWedges << std::endl;
	} else {
		std::cerr << "sim histogram for rings and/or wedges has zero integral and thus cannot be scaled" << std::endl;
		return;
	}


	const int rebinfactor = 4;

	double maxdatarings = hDataRings->GetMaximum();
	double maxdatawedges = hDataWedges->GetMaximum();

	double maxsimrings = hSimRings->GetMaximum();
	double maxsimwedges = hSimWedges->GetMaximum();

	double ymaxrings = std::max(maxdatarings,maxsimrings)*1.1*rebinfactor;
	double ymaxwedges = std::max(maxdatawedges,maxsimwedges)*1.1*rebinfactor;

	hDataRings->SetMaximum(ymaxrings);
	hDataWedges->SetMaximum(ymaxwedges);

	//start drawing w/ canvas

	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);

	//calculate chi squareds here:
	double chi2_rings = hDataRings->Chi2Test(hSimRings, "CHI2/NDF UW");
	double chi2_wedges = hDataWedges->Chi2Test(hSimWedges, "CHI2/NDF UW");


	hDataRings->SetLineColor(kViolet);
	hDataRings->SetLineWidth(4);
	TString title = "^{6}Li (0^{+}) Data vs Sim (SABRE3 Ring Summary);Energy (MeV)";
	hDataRings->SetTitle(title);
	TH1D* hDataRings_rebin = dynamic_cast<TH1D*>(hDataRings->Rebin(rebinfactor));
	hDataRings_rebin->GetXaxis()->SetRangeUser(0,4);
	hDataRings_rebin->GetYaxis()->SetRangeUser(0,200);
	hDataRings_rebin->Draw("HIST");

	hSimRings->SetLineColor(kOrange);
	hSimRings->SetLineWidth(2);
	TH1D* hSimRings_rebin = dynamic_cast<TH1D*>(hSimRings->Rebin(rebinfactor));
	hSimRings_rebin->GetXaxis()->SetRangeUser(0,4);
	hSimRings_rebin->GetYaxis()->SetRangeUser(0,200);
	hSimRings_rebin->Draw("HIST SAME");

	//add legend
	TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hDataRings,"Data","l");
	TString simringlabel = Form("Sim (scaled x %.3f)",scaleFactorRings);
	legend->AddEntry(hSimRings,simringlabel,"l");
	legend->Draw();

	//add text box to differentiate pngs!
	TPaveText *pt = new TPaveText(0.65,0.63,0.88,0.74,"NDC");
	pt->AddText(beamstring);
	pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2_rings));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();


	c1->Update();
	//save to png file:
	c1->SaveAs(Form("Lithium6_0plus_simVSdata_SABRE3_RingSummary_%s.png",beamstring.Data()));


	//clear canvas, do same thing but for wedges this time!
	c1->Clear();

	hDataWedges->SetLineColor(kViolet);
	hDataWedges->SetLineWidth(4);
	title = "^{6}Li (0^{+}) Data vs Sim (SABRE3 Wedge Summary);Energy (MeV)";
	hDataWedges->SetTitle(title);
	TH1D* hDataWedges_rebin = dynamic_cast<TH1D*>(hDataWedges->Rebin(rebinfactor));
	hDataWedges_rebin->GetXaxis()->SetRangeUser(0,4);
	hDataWedges_rebin->GetYaxis()->SetRangeUser(0,200);
	hDataWedges_rebin->Draw("HIST");

	hSimWedges->SetLineColor(kOrange);
	hSimWedges->SetLineWidth(2);
	TH1D* hSimWedges_rebin = dynamic_cast<TH1D*>(hSimWedges->Rebin(rebinfactor));
	hSimWedges_rebin->GetXaxis()->SetRangeUser(0,4);
	hSimWedges_rebin->GetYaxis()->SetRangeUser(0,200);
	hSimWedges_rebin->Draw("HIST SAME");

	//add legend:
	legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hDataWedges,"Data","l");
	simringlabel = Form("Sim (scaled x %.3f)",scaleFactorWedges);
	legend->AddEntry(hSimWedges,simringlabel,"l");
	legend->Draw();

	//add text box to differentiate pngs!
	pt = new TPaveText(0.65,0.63,0.88,0.74,"NDC");
	pt->AddText(beamstring);
	pt->AddText(Form("#chi^{2}/NDF = %.3f",chi2_wedges));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	c1->Update();
	//save to png file:
	c1->SaveAs(Form("Lithium6_0plus_simVSdata_SABRE3_WedgeSummary_%s.png",beamstring.Data()));

}

void Lithium6_0plus_sabre3summaries_auto(){
	std::vector<TString> beamstrings = {
										"fixedpoint",
										"gaus001",
										"gaus002",
										"gaus003",
										"gaus004",
										"gaus005",
										"gaus006",
										"gaus007",
										"gaus008",
										"gaus009",
										"gaus010"
									};

	for(const auto& bs : beamstrings){
		Lithium6_0plus_sabre3summaries(bs);
	}
}

void Lithium6_3plus_1par(){

}

void Lithium6_3plus_2par(){
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/NEW_LiFha_2par_exp_3plus_OUTPUT_FORCETHETA20PHI0.root";
	TString dataHistLocalPath = "2par/3plus/hSABRE_ChannelHits_3plus";

	TString simFilePath = "/mnt/e/SABREsim/det/kin3mc/kin3mc_7Li3He4He6Li2186keV_4He2H_7500keV_moretheta_offset0.root";
	TString simHistLocalPath = "2par/3plus/h2SABRE_ChannelHits_3plus";

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

	if(hData->GetNbinsX() != hSim->GetNbinsX()){
		std::cerr << "Histogram bin counts do not match!" << std::endl;
		return;
	}

	double dataIntegral = hData->Integral();
	double simIntegral = hSim->Integral();

	double scaleFactor = 0.;
	if(simIntegral > 0){
		scaleFactor = dataIntegral/simIntegral;
		hSim->Scale(scaleFactor);
		std::cout << "Applied global scale factor " << scaleFactor << " to sim" << std::endl;
	} else {
		std::cerr << "sim histogram has zero integral and thus cannot scale" << std::endl;
		return;
	}

	TCanvas *c1 = new TCanvas("c1", "Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);
	const int rebinfactor = 1;

	double maxdata = hData->GetMaximum();
	double maxsim = hSim->GetMaximum();

	double ymax = std::max(maxdata,maxsim)*1.1;
	hData->SetMaximum(ymax);

	hData->SetLineColor(kViolet);
	hData->SetLineWidth(4);
	hData->SetTitle("^{6}Li (3^{+} E = 2.186 MeV) Data vs Sim;Energy (MeV);Counts");
	hData->Draw("HIST");

	hSim->SetLineColor(kOrange);
	hSim->SetLineWidth(4);
	hSim->Draw("HIST SAME");

	TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hData,"Data","l");
	TString simlabel = Form("Sim (scaled x %.3f)", scaleFactor);
	legend->AddEntry(hSim,simlabel,"l");
	legend->Draw();

	c1->Update();

	//save to a new root file:
	TString outfile_root_name = "Lithium6_3plus_2par_simVSdata.root";
	TFile *outfile_root = new TFile(outfile_root_name,"RECREATE");
	if(!outfile_root || outfile_root->IsZombie()){
		std::cerr << "Error creating output root file" << std::endl;
		return;
	}

	hData->Write("hData");
	hSim->Write("hSim");
	c1->Write("overlay_canvas");
	outfile_root->Close();

	std::cout << "Finished! Output written to: " << outfile_root_name << std::endl;

}


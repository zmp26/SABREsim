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
#include <TH2I.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TPaveText.h>
#include <cmath>
#include <TMath.h>
#include <string>

bool isInt(double x, double epsilon=1e-6){
	return std::fabs(x - std::round(x)) < epsilon;
}

int neededPrecision(double x, double epsilon=1e-6){
	double frac = std::fabs(x - std::round(x));
	// if(frac == 0.0) return 1;
	// if(std::fmod(frac*10,1.0) == 0.) return 1;
	// if(std::fmod(frac*100,1.0) == 0.) return 2;
	// if(std::fmod(frac*1000,1.0) == 0.) return 3;
	// return 3;

	if(isInt(frac)) return 1;
	if(isInt(frac*10)) return 1;
	if(isInt(frac*100)) return 2;
	if(isInt(frac*1000)) return 3;
	return 3;
}

std::string smartFormat(double x){
	int prec = neededPrecision(x);
	std::ostringstream ss;
	ss << std::fixed << std::setprecision(prec) << x;
	return ss.str();
}

double ComputeChi2BinByBin(TH1* hData, TH1* hSim, double &ndf){

	if(!hData || !hSim){
		std::cerr << "Error: null histo pointer for hData or hSim" << std::endl;
		ndf = 0;
		return -1;
	}

	if(hData->GetNbinsX() != hSim->GetNbinsX()){
		std::cerr << "Histograms have inequal binning!" << std::endl;
		ndf = 0;
		return -1;
	}

	double chi2 = 0.;
	ndf = 0.;

	int nbins = hData->GetNbinsX();
	for(int i=1; i<=nbins; i++){
		double dataval = hData->GetBinContent(i);
		double simval = hSim->GetBinContent(i);

		double datasigma2 = (dataval > 0) ? dataval : 1.;//avoids dividing by 0 error for 0 count bins -> assumes poisson statistics
		double simsigma2 = (simval > 0) ? simval : 1.;//avoids dividing by 0 error, assumes poisson stats

		double numerator = std::pow((dataval - simval), 2);
		double denominator = datasigma2 + simsigma2;
		double tempchi2 = (numerator/denominator);

		if(dataval != 0 || simval != 0){
			ndf += 1.;
			chi2 += tempchi2;
		}

	}

	return chi2;

}

double ComputeChi2BinByBin(TH1* hData, TH1* hSim_unscaled, double scaleFactor, double &ndf){
	//same function as ComputeChi2BinByBin but allows for passing the unscaled sim histogram and thus the better relative sigmas when scaled!

	if(!hData || !hSim_unscaled){
		std::cerr << "Error: null histo pointer for hData or hSim_unscaled" << std::endl;
		ndf = 0;
		return -1;
	}

	if(hData->GetNbinsX() != hSim_unscaled->GetNbinsX()){
		std::cerr << "Histograms have inequal binning!" << std::endl;
		ndf = 0;
		return -1;
	}

	double chi2 = 0.;
	ndf = 0.;

	//scale hSim_unscaled by scalefactor:
	hSim_unscaled->Scale(scaleFactor);

	//now iterate through bins:
	int nbins = hData->GetNbinsX();
	for(int i = 1; i <= nbins; i++){
		double dataval = hData->GetBinContent(i);
		double dataerror = hData->GetBinError(i);//sigma

		double simval = hSim_unscaled->GetBinContent(i);
		double simerror = hSim_unscaled->GetBinError(i);//sigma

		double sigma2 = dataerror*dataerror + simerror*simerror;//add in quadrature!

		double numerator = std::pow((dataval - simval), 2);
		double denominator = sigma2;

		double tempchi2 = (numerator/denominator);

		if(dataval != 0 || simval != 0){
			ndf += 1.;
			chi2 += tempchi2;
		} 

	}

	return chi2;

}

void Lithium6_1plus(int ring, TString beamstringx, TString beamstringy){

	TString anglestring = "16782278";

	//uncomment below for surface laptop
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = Form("1par/1plus/hSABRE_SABRE3_Ring%dESummary_1plus",ring);

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%sx_%sy_histos.root", anglestring.Data(), beamstringx.Data(), beamstringy.Data());
	TString simHistLocalPath = Form("SABRE/SABRE3/Summary/hSABRE3_Ring%dSummary",ring);

	//std::cout << simFilePath << std::endl;

	//uncomment below for DESKTOP
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = Form("1par/1plus/hSABRE_SABRE3_Ring%dESummary_1plus",ring);

	// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%sx_%sy_histos.root", anglestring.Data(), beamstringx.Data(), beamstringy.Data());
	// TString simHistLocalPath = Form("SABRE/SABRE3/Summary/hSABRE3_Ring%dSummary",ring);

	//for our sps acceptances, rings 7 and 8 on SABRE3 illuminate
	// if(ring != 7 && ring != 8) {
	// 	std::cout << "ring not found, try 7 or 8" << std::endl;
	// 	return;
	// }

	//prep output file name using ring:
	TString outfile_root_name;
	//outfile_root_name = Form("Lithium6_1plus_simVSdata_ring%d.root",ring);//"Lithium6_0plus_simVSdata_ELoss2.root";
	outfile_root_name = Form("Lithium6_1plus_simVSdata_%sx_%sy_theta%s_ring%d.root", beamstringx.Data(), beamstringy.Data(), anglestring.Data(), ring);



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
	pt->AddText(Form("%sx %sy SPSTheta%s", beamstringx.Data(), beamstringy.Data(), anglestring.Data()));
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
	c1->SaveAs(Form("Lithium6_1plus_simVSdata_%sx_%sy_theta%s_ring%d.png", beamstringx.Data(), beamstringy.Data(), anglestring.Data(), ring));
	//outfile_root->Close();

	//std::cout << "Finished! Output written to: " << outfile_root_name << std::endl;
	//std::cout << "Scale factors saved to Lithium6_1plus_scales.txt" << std::endl;
}

void Lithium6_1plus_auto(){
	std::vector<int> rings = {7,8,9};
	std::vector<TString> beamstringsx = {
										// "fixedpoint",
										// "gaus001",
										// "gaus002",
										// "gaus003",
										// "gaus004",
										// "gaus005"
										// "gaus006",
										// "gaus007",
										// "gaus008",
										// "gaus009",
										// "gaus010"
										"gaus004",
									};
	std::vector<TString> beamstringsy = {"gaus0015"};

	for(size_t r=0; r<rings.size(); r++){
		for(size_t bx=0; bx<beamstringsx.size(); bx++){
			//Lithium6_1plus(rings[r],beamstrings[b]);
			for(size_t by=0; by<beamstringsy.size(); by++){
				Lithium6_1plus(rings[r], beamstringsx[bx], beamstringsy[by]);
			}
		}
	}
}

TString Lithium6_1plus_sabrehits(TString beamstringx, TString beamstringy, TString anglestring, TString phistring, TString stragglestring){

	//TString anglestring = "16782278";

	//uncomment below for DESKTOP
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%sx_%sy_histos.root",anglestring.Data(),beamstringx.Data(),beamstringy.Data());
	// TString simHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	//uncomment below for LAPTOP
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_phi_%s_%sx_%sy_%s_histos.root",anglestring.Data(), phistring.Data(), beamstringx.Data(),beamstringy.Data(),stragglestring.Data());
	TString simHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	TString pngoutname = Form("Lithium6_1plus_simVSdata_sabrehits_theta%s_phi_%s_%sx_%sy_%s.png",anglestring.Data(), phistring.Data(), beamstringx.Data(), beamstringy.Data(), stragglestring.Data());


	//prep output file name
	//TString outfile_root_name = Form("Lithium6_1plus_simVSdata_sabrehits_%s.root",beamstring.Data());


	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return "ERROR";
	}

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram" << std::endl;
		datafile->Close();
		return "ERROR";
	}
	hData->SetDirectory(0);
	datafile->Close();



	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return "ERROR";
	}

	TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
	if(!hSim){
		std::cerr << "Error retrieving sim histogram" << std::endl;
		simfile->Close();
		return "ERROR";
	}
	hSim->SetDirectory(0);
	simfile->Close();



	if((hData->GetNbinsX() != hSim->GetNbinsX()) || (hData->GetNbinsY() != hSim->GetNbinsY())){
		std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
		return "ERROR";
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
		return "ERROR";
	}

	//calculate chi2:
	double chi2 = hData->Chi2Test(hSim, "CHI2/NDF UW");


	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);



	hSim->Draw();

	TPaveText *pt = new TPaveText(0.65, 0.75, 0.88, 0.88, "NDC");
	pt->AddText(Form("%sx %sy (scale factor = %.3f)", beamstringx.Data(), beamstringy.Data(), scaleFactor));
	pt->AddText(Form("#chi^{2}/NDF = %.1f",chi2));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	TPaveText *ptangle = new TPaveText(0.12, 0.7, 0.27, 0.9, "NDC");
	TString temp;
	if(anglestring == "188208") temp = "19.8#circ #pm 1#circ";
	else if(anglestring == "178218") temp = "19.8#circ #pm 2#circ";
	else if(anglestring == "168228") temp = "19.8#circ #pm 3#circ";
	else if(anglestring == "158238") temp = "19.8#circ #pm 4#circ";
	else if(anglestring == "148248") temp = "19.8#circ #pm 5#circ";
	ptangle->AddText(Form("%s",temp.Data()));
	ptangle->AddText(Form("phi %s", phistring.Data()));
	ptangle->AddText(stragglestring);
	//ptangle->AddText(Form("#chi^{2}/NDF = %.1f",chi2));
	ptangle->SetFillColorAlpha(kWhite,0.5);
	ptangle->Draw();

	c1->Update();

	// TFile* outfile_root = new TFile(outfile_root_name,"RECREATE");
	// if(!outfile_root || outfile_root->IsZombie()){
	// 	std::cerr << "Error creating output file" << std::endl;
	// 	return;
	// }

	//hData->Write("hData_original");
	//hSim->Write("hSim_scaled");
	c1->SaveAs(pngoutname);
	c1->Clear();
	hData->Draw();
	pt = new TPaveText(0.65, 0.75, 0.88, 0.88, "NDC");
	pt->AddText("LiFha exp 1par ^{6}Li 1^{+} DATA");
	c1->Update();
	c1->SaveAs(Form("Lithium6_1plus_simVSdata_sabrehits_%s.png","1PAREXPDATA"));
	//outfile_root->Close();
	return pngoutname;
	//std::cout << "Finished! Output written to: " << outfile_root_name << std::endl;

}

void Lithium6_1plus_sabrehits_auto(){

	// std::vector<TString> beamstrings = {
	// 									"fixed",
	// 									"gaus001",
	// 									"gaus002",
	// 									"gaus003",
	// 									"gaus004",
	// 									"gaus005"
	// 									// "gaus006",
	// 									// "gaus007",
	// 									// "gaus008",
	// 									// "gaus009",
	// 									// "gaus010"
	// 								};

	std::vector<TString> beamstringsx = {
											"gaus001"
										};

	std::vector<TString> beamstringsy = {
											"gaus001"									
										};

	std::vector<TString> anglestrings = {
											//"193203",			// 19.8 +/- 0.5
											//"188208",			// 19.8 +/- 1.0
											"178218"			// 19.8 +/- 2.0
											//"168228",			// 19.8 +/- 3.0
											//"158238",			// 19.8 +/- 4.0
											//"148248"			// 19.8 +/- 5.0
										};

	TString straggle = "straggleOn";

	//std::vector<double> phis = {0.5, 1.0, 1.5, 2.0, 2.5};
	std::vector<double> phis = {1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0};
	std::vector<TString> phistrings;
	for(const auto& philower : phis){
		for(const auto& phiupper : phis){
			//phistrings.push_back(Form("%f_%f",-philower,phiupper));
			TString s = TString::Format("%s_%s", smartFormat(-philower).c_str(), smartFormat(phiupper).c_str());
			phistrings.push_back(s);
		}
	}



	for(const auto& bsx : beamstringsx){
		for(const auto& bsy : beamstringsy){
			std::vector<TString> pngs;
			for(const auto& angstr : anglestrings){
				for(const auto& phistr : phistrings){
					TString pngoutname = Lithium6_1plus_sabrehits(bsx, bsy, angstr, phistr, straggle);
					pngs.push_back(pngoutname);
				}
			}

			//for all anglestring sabre hit patterns for this given bsx, bsy combo let's make into a gif:
			TString gifname = Form("Lithium6_1plus_simVSdata_sabrehits_%sx_%sy.gif",bsx.Data(),bsy.Data());
			TString cmd = "convert -delay 100 -loop 0 ";
			for(const auto& fn : pngs){
				cmd += fn;
				cmd += " ";
			}
			cmd += gifname;

			std::cout << "Creating gif for bsx = " << bsx << " and bsy = " << bsy << std::endl;
			gSystem->Exec(cmd);
			std::cout << "Saved " << gifname << std::endl;

			//clean up pngs:
			cmd = "rm ";
			for(const auto& fn : pngs){
				cmd += fn;
				cmd += " ";
			}

			// gSystem->Exec(cmd);
		}
	}

}

void Lithium6_1plus_sabre3summaries(TString beamstring){

	TString anglestring = "16782278";

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

TString Lithium6_1plus_pixelhistos(int ringChan, int wedgeChan, TString beamstringx, TString beamstringy, TString anglestring, TString phistring, TString stragglestring){

	int SABRE_ID = 3;

	TString pixelhistoname = Form("hSABRE%d_pixel_r%dw%d_ESummary",SABRE_ID,ringChan,wedgeChan);
	//TString anglestring = "16782278";

	//uncomment below for DESKTOP:
	//TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	//TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%s_histos.root",anglestring.Data(),beamstring.Data());

	//uncomment below for LAPTOP:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = Form("SABRE/SABRE%d/Pixels/",SABRE_ID);
	dataHistLocalPath = dataHistLocalPath + pixelhistoname;

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_phi_%s_%sx_%sy_%s_histos.root",anglestring.Data(), phistring.Data(), beamstringx.Data(), beamstringy.Data(), stragglestring.Data());
	TString simHistLocalPath = Form("SABRE/SABRE%d/Pixels/",SABRE_ID);
	simHistLocalPath = simHistLocalPath + pixelhistoname;


	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		std::cerr << "dataFilePath = " << dataFilePath << "\n\n";
		return "ERROR";
	}

	TH1 *hData = dynamic_cast<TH1*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram" << std::endl;
		std::cerr << "dataFilePath = " << dataFilePath << "\n";
		std::cerr << "dataHistLocalPath = " << dataHistLocalPath << "\n\n";
		datafile->Close();
		return "ERROR";
	}
	hData->SetDirectory(0);
	datafile->Close();



	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return "ERROR";
	}

	TH1 *hSim = dynamic_cast<TH1*>(simfile->Get(simHistLocalPath));
	if(!hSim){
		std::cerr << "Error retrieving sim histogram" << std::endl;
		simfile->Close();
		return "ERROR";
	}
	hSim->SetDirectory(0);
	simfile->Close();



	if((hData->GetNbinsX() != hSim->GetNbinsX()) || (hData->GetNbinsY() != hSim->GetNbinsY())){
		std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
		return "ERROR";
	}


	double dataIntegral = hData->Integral();
	double simIntegral = hSim->Integral();

	double scaleFactor = 0.;

	if(simIntegral > 0){
		scaleFactor = dataIntegral/simIntegral;
		//hSim->Scale(scaleFactor);
		std::cout << "Applied global scale factor to sim: " << scaleFactor << std::endl;
	} else {
		std::cerr << "sim histogram has zero integral and thus cannot scale" << std::endl;
		return "ERROR";
	}

	const int rebinfactor = 4;

	double maxdata = hData->GetMaximum();
	double maxsim = hSim->GetMaximum();

	double ymax = std::max(maxdata,maxsim);

	//calculate chi2:
	//double chi2 = hData->Chi2Test(hSim, "CHI2/NDF UW");
	double ndf;
	double chi2 = ComputeChi2BinByBin(hData, hSim, scaleFactor, ndf);//ComputeChi2BinByBin(TH1* hData, TH1* hSim_unscaled, double scaleFactor, double &ndf)


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
	//pt->AddText(Form("%s",anglestring.Data()));
	pt->AddText(Form("%sx %sy", beamstringx.Data(), beamstringy.Data()));
	pt->AddText(Form("(scale factor = %.3f)", scaleFactor));
	pt->AddText(Form("#chi^{2}/NDF = %.1f",chi2/ndf));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	TPaveText *ptangle = new TPaveText(0.12, 0.7, 0.27, 0.9, "NDC");
	TString temp;
	if(anglestring == "193203") temp = "19.8#circ #pm 0.5#circ";
	else if(anglestring == "188208") temp = "19.8#circ #pm 1#circ";
	else if(anglestring == "178218") temp = "19.8#circ #pm 2#circ";
	else if(anglestring == "168228") temp = "19.8#circ #pm 3#circ";
	else if(anglestring == "158238") temp = "19.8#circ #pm 4#circ";
	else if(anglestring == "148248") temp = "19.8#circ #pm 5#circ";
	ptangle->AddText(Form("%s",temp.Data()));
	ptangle->AddText(stragglestring);
	//ptangle->AddText(Form("#chi^{2}/NDF = %.1f",chi2));
	ptangle->SetFillColorAlpha(kWhite,0.5);
	ptangle->Draw();

	c1->Update();

	TString pngoutname = Form("Lithium6_1plus_simVSdata_pixelhisto_r%dw%d_theta%s_%sx_%sy_%s.png",ringChan,wedgeChan,anglestring.Data(),beamstringx.Data(), beamstringy.Data(), stragglestring.Data());
	c1->SaveAs(pngoutname);

	delete c1;
	delete hData;
	//delete hDataRebin;
	delete hSim;
	//delete hSimRebin;

	return pngoutname;
}

void Lithium6_1plus_pixelhistos_auto(){

	TString beamx = "gaus0005";
	TString beamy = "gaus0005";

	TString anglestring = "178218";

	TString phistring = "-2.125_2.125";

	TString straggle = "straggleOff";

	std::vector<int> ringChans = {
									71,
									72,
									73	
								 };

	std::vector<int> wedgeChans = {
									29,
									30
								  };

	for(const auto& ring : ringChans){
		for(const auto& wedge : wedgeChans){

			Lithium6_1plus_pixelhistos(ring, wedge, beamx, beamy, anglestring, phistring, straggle);

		}
	}

	std::cout << "Done!" << std::endl;

}

// void Lithium6_1plus_pixelhistos_auto(){
// 	// std::vector<TString> beamstrings = {
// 	// 									"fixedpoint",
// 	// 									"gaus001",
// 	// 									"gaus002",
// 	// 									"gaus003",
// 	// 									"gaus004",
// 	// 									"gaus005",
// 	// 									// "gaus006",
// 	// 									// "gaus007",
// 	// 									// "gaus008",
// 	// 									// "gaus009",
// 	// 									// "gaus010"
// 	// 								};

// 	TString beamx = "gaus0005";
// 	TString beamy = "gaus0005";

// 	std::vector<TString> anglestrings = {"178218"};

// 	// //std::vector<TString> anglestrings = {
// 	// 										"193203",			// 19.8 +/- 0.5
// 	// 										"188208",			// 19.8 +/- 1.0
// 	// 										"178218",			// 19.8 +/- 2.0
// 	// 										"168228",			// 19.8 +/- 3.0
// 	// 										"158238",			// 19.8 +/- 4.0
// 	// 										"148248"			// 19.8 +/- 5.0
// 	// 									};

// 	// x: fixed 			gaus0005
// 	// y: fixed 			gaus00025

// 	std::vector<int> ringChans = {
// 									71,
// 									72,
// 									73	
// 								 };

// 	std::vector<int> wedgeChans = {
// 									29,
// 									30
// 								  };

// 	// for(const auto& bs : beamstrings){
// 	// 	for(const auto& ring : ringChans){
// 	// 		for(const auto& wedge : wedgeChans){
// 	// 			Lithium6_1plus_pixelhistos(ring,wedge,bs);
// 	// 		}
// 	// 	}
// 	// }

// 	for(const auto& ring : ringChans){
// 		for(const auto& wedge : wedgeChans){

// 			std::vector<TString> pngs;
			
// 			for(const auto& angstr : anglestrings){
// 				TString pngoutname = Lithium6_1plus_pixelhistos(ring,wedge,beamx,beamy,angstr);
// 				pngs.push_back(pngoutname);
// 			}

// 			//we have now run Lithium6_1plus_pixelhistos for all angstrs for this (ring,wedge) pair, so combine into gif!

// 			//combing all (r,w) pair pngs into gif using system call to image magick:
// 			TString gifname = Form("r%dw%d_%sx_%sy.gif",ring,wedge,beamx.Data(),beamy.Data());
// 			TString cmd = "convert -delay 100 -loop 0 ";
// 			for(const auto& fn : pngs){
// 				cmd += fn;
// 				cmd += " ";
// 			}
// 			cmd += gifname;

// 			std::cout << "Creating gif for (r,w) = (" << ring << ", " << wedge << ")\n";
// 			//std::cout << "\tcmd = " << cmd << "\n\n";
// 			gSystem->Exec(cmd);
// 			std::cout << "Saved " << gifname << std::endl;

// 			//clean up pngs:
// 			cmd = "rm ";
// 			for(const auto& fn : pngs){
// 				cmd += fn;
// 				cmd += " ";
// 			}

// 			gSystem->Exec(cmd);
// 		}
// 	}

// }

std::pair<TString, TString> parsephistring(const TString& ts){
	std::string s = ts.Data();
	auto pos = s.find('_');
	if(pos == std::string::npos) throw std::runtime_error("Invalid phi string: " + s);

	TString low = s.substr(0,pos);
	TString high = s.substr(pos+1);

	return {low, high};
}

void Lithium6_1plus_fourpixelchisquared(){

	TString anglestring = "178218";

	std::vector<TString> beamstrings = {
											//"fixed"
											"gaus001"
											// "gaus002",
											// "gaus003",
											// "gaus004",
											// "gaus005"
										};

	// TString bsx = "fixed";
	// TString bsy = "fixed";


	//std::vector<double> phis = {1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0};
	std::vector<double> phis = {2.0, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3.0, 3.125, 3.25, 3.375, 3.5};
	std::vector<TString> phistrings;
	for(const auto& philower : phis){
		for(const auto& phiupper : phis){
			//phistrings.push_back(Form("%f_%f",-philower,phiupper));
			TString s = TString::Format("%s_%s", smartFormat(-philower).c_str(), smartFormat(phiupper).c_str());
			phistrings.push_back(s);
		}
	}

	double xedges[14] = {-3.5625, -3.4375, -3.3125, -3.1875, -3.0625, -2.9375, -2.8125, -2.6875, -2.5625, -2.4375, -2.3125, -2.1875, -2.0625, -1.9375};
	double yedges[14] = {1.9375, 2.0625, 2.1875, 2.3125, 2.4375, 2.5625, 2.6875, 2.8125, 2.9375, 3.0625, 3.1875, 3.3125, 3.4375, 3.5625};

	TH2D *hGridSearchChi2 = new TH2D("hGridSearchChi2", "GridSearchChi2", 13, xedges, 13, yedges);
	TH2D *hGridSearchReducedChi2 = new TH2D("hGridSearchReducedChi2", "GridSearchReducedChi2", 13, xedges, 13, yedges);

	//---------------------------------------------------
	//				Establish data values
	//---------------------------------------------------

	//uncomment for DESKTOP:
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
	//uncomment for LAPTOP:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";


	// TString path_pix_r71_w29 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r71_w29";
	// TString path_pix_r72_w29 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r72_w29";
	// TString path_pix_r72_w30 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r72_w30";
	// TString path_pix_r71_w30 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r71_w30";

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	// TH1 *hpix_r71_w29 = dynamic_cast<TH1*>(datafile->Get(path_pix_r71_w29));
	// TH1 *hpix_r72_w29 = dynamic_cast<TH1*>(datafile->Get(path_pix_r72_w29));
	// TH1 *hpix_r72_w30 = dynamic_cast<TH1*>(datafile->Get(path_pix_r72_w30));
	// TH1 *hpix_r71_w30 = dynamic_cast<TH1*>(datafile->Get(path_pix_r71_w30));
	// if(!hpix_r71_w29 || !hpix_r72_w29 || !hpix_r72_w30 || !hpix_r71_w30){
	// 	std::cerr << "Error retrieving at least one data histogram!" << std::endl;
	// 	return;
	// }

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram!" << std::endl;
		return;
	}

	// hpix_r71_w29->SetDirectory(0);
	// hpix_r72_w29->SetDirectory(0);
	// hpix_r72_w30->SetDirectory(0);
	// hpix_r71_w30->SetDirectory(0);
	hData->SetDirectory(0);
	datafile->Close();

	// double counts_pix_r71_w29 = hpix_r71_w29->GetEntries();
	// double counts_pix_r72_w29 = hpix_r72_w29->GetEntries();
	// double counts_pix_r72_w30 = hpix_r72_w30->GetEntries();
	// double counts_pix_r71_w30 = hpix_r71_w30->GetEntries();

	double dataIntegral = (hData->GetBinContent(hData->GetBin(30,8))) + (hData->GetBinContent(hData->GetBin(30,9))) + (hData->GetBinContent(hData->GetBin(31,9))) + (hData->GetBinContent(hData->GetBin(31,8)));

	//std::vector<double> fourpix_relcounts = {counts_pix_r71_w29/fourpixsum, counts_pix_r72_w29/fourpixsum, counts_pix_r72_w30/fourpixsum, counts_pix_r71_w30/fourpixsum};

	//for(const auto& fn : filenames){

	for(const auto& ps : phistrings){
		for(const auto& bsx : beamstrings){
			for(const auto& bsy : beamstrings){

				TString fn = Form("kin2mc_7Li3He4He6Ligs_7500keV_theta%s_phi_%s_%sx_%sy_histos.root", anglestring.Data(), ps.Data(), bsx.Data(), bsy.Data());

				//establish sim values:
				//uncomment for DESKTOP
				// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/%s",fn.Data());
				// TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
				//uncomment for LAPTOP:
				TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/%s",fn.Data());
				TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";


				TFile *simfile = new TFile(simFilePath,"READ");
				if(!simfile || simfile->IsZombie()){
					std::cerr << "Error opening sim file " << simFilePath << "\n\n";
					continue;
				}


				TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
				if(!hSim){
					std::cerr << "Error retrieving sim histo from " << simFilePath << "\n\n";
					continue;
				}

				double simIntegral = (hSim->GetBinContent(hSim->GetBin(30,8))) + (hSim->GetBinContent(hSim->GetBin(30,9))) + (hSim->GetBinContent(hSim->GetBin(31,9))) + (hSim->GetBinContent(hSim->GetBin(31,8)));


				hSim->SetDirectory(0);
				simfile->Close();

				if( hSim->GetNbinsX() != hData->GetNbinsX() ){
					std::cerr << "Histogram binning does not match for data/sim histogram pair\n";
					continue;
				}




				double scaleFactor;
				if(simIntegral > 0){
					scaleFactor = dataIntegral/simIntegral;
					hSim->Scale(scaleFactor);
				} else {
					std::cerr << "sim histogram from " << fn.Data() << " has zero integral and thus cannot be scaled! Continuing..." << std::endl;
					continue;
				}

				//calculate the chi squared for this sim file:

				double r71w29_data = (hData->GetBinContent(hData->GetBin(30,8)));//binx = 30 is wedge 29, biny = 8 is local ring 7 which for SABRE3 is ring 71
				double r72w29_data = (hData->GetBinContent(hData->GetBin(30,9)));//binx = 30 is wedge 29, biny = 9 is local ring 8 which for SABRE3 is ring 72
				double r72w30_data = (hData->GetBinContent(hData->GetBin(31,9)));//binx = 31 is wedge 30, biny = 9 is local ring 8 which for SABRE3 is ring 72
				double r71w30_data = (hData->GetBinContent(hData->GetBin(31,8)));//binx = 31 is wedge 30, biny = 8 is local ring 7 which for SABRE3 is ring 71

				double r71w29_data_error = (hData->GetBinError(hData->GetBin(30,8)));
				double r72w29_data_error = (hData->GetBinError(hData->GetBin(30,9)));
				double r72w30_data_error = (hData->GetBinError(hData->GetBin(31,9)));
				double r71w30_data_error = (hData->GetBinError(hData->GetBin(31,8)));

				double r71w29_sim = (hSim->GetBinContent(hSim->GetBin(30,8)));
				double r72w29_sim = (hSim->GetBinContent(hSim->GetBin(30,9)));
				double r72w30_sim = (hSim->GetBinContent(hSim->GetBin(31,9)));
				double r71w30_sim = (hSim->GetBinContent(hSim->GetBin(31,8)));

				double r71w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,8)));
				double r72w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,9)));
				double r72w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,9)));
				double r71w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,8)));

				double chi2 = ( std::pow((r71w29_data - r71w29_sim),2) / ( std::pow(r71w29_data_error,2) + std::pow(r71w29_sim_error,2)) ) + 
							  ( std::pow((r72w29_data - r72w29_sim),2) / ( std::pow(r72w29_data_error,2) + std::pow(r72w29_sim_error,2)) ) +
							  ( std::pow((r72w30_data - r72w30_sim),2) / ( std::pow(r72w30_data_error,2) + std::pow(r72w30_sim_error,2)) ) +
							  ( std::pow((r71w30_data - r71w30_sim),2) / ( std::pow(r71w30_data_error,2) + std::pow(r71w30_sim_error,2)) );

				Int_t ndf = 4 - 1;//4 pixels minus the 1 DOF (the scale factor is only fitted parameter here)

				double reducedchi2 = chi2/ndf;


				//fill chi2 histogram from grid search:
				int xbin = 0;
				int ybin = 0;

				// if(bsx == "fixed") xbin = 0;
				// else if(bsx == "gaus001") xbin = 1;
				// else if(bsx == "gaus002") xbin = 2;
				// else if(bsx == "gaus003") xbin = 3;
				// else if(bsx == "gaus004") xbin = 4;
				// else if(bsx == "gaus005") xbin = 5;

				// if(bsy == "fixed") ybin = 0;
				// else if(bsy == "gaus001") ybin = 1;
				// else if(bsy == "gaus002") ybin = 2;
				// else if(bsy == "gaus003") ybin = 3;
				// else if(bsy == "gaus004") ybin = 4;
				// else if(bsy == "gaus005") ybin = 5;

				std::pair<TString,TString> philowhigh = parsephistring(ps);
				TString philow = philowhigh.first;
				TString phihigh = philowhigh.second;

				// if(philow == "-2.5") xbin = 5;
				// else if(philow == "-2.0") xbin = 4;
				// else if(philow == "-1.5") xbin = 3;
				// else if(philow == "-1.0") xbin = 2;
				// else if(philow == "-0.5") xbin = 1;

				// if(phihigh == "2.5") ybin = 5;
				// else if(phihigh == "2.0") ybin = 4;
				// else if(phihigh == "1.5") ybin = 3;
				// else if(phihigh == "1.0") ybin = 2;
				// else if(phihigh == "0.5") ybin = 1;

				hGridSearchChi2->Fill(philow.Atof(),phihigh.Atof(),chi2);
				hGridSearchReducedChi2->Fill(philow.Atof(),phihigh.Atof(),reducedchi2);
			}
		}
	}

	//add labels showing value and relative rank on each bin:
	int nx = hGridSearchReducedChi2->GetNbinsX();
	int ny = hGridSearchReducedChi2->GetNbinsY();

	//store (value, global_bin) for ranking
	std::vector<std::pair<double,int>> values;
	values.reserve(nx*ny);

	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){
			double val = hGridSearchReducedChi2->GetBinContent(ix,iy);
			int globalbin = hGridSearchReducedChi2->GetBin(ix,iy);

			if(val > 0) values.push_back({val,globalbin});
		}
	}

	//sort by reduced chi2 value:
	std::sort(values.begin(), values.end(),
		      [](auto&a, auto&b){ return a.first < b.first; });

	std::map<int,int> rankmap;
	for(size_t i=0; i<values.size(); i++){
		rankmap[values[i].second] = i+1;
	}

	TLatex latex;
	latex.SetTextAlign(22);
	latex.SetTextSize(0.012);
	latex.SetTextColor(kBlack);

	TFile *outfile = new TFile("GridSearch.root","RECREATE");

	outfile->cd();
	//hData->Write();
	//hSim->Write();
	hGridSearchChi2->Write();
	hGridSearchReducedChi2->Write();

	//set up square TCanvas here:
	TCanvas *c1 = new TCanvas("c1","Grid Search Reduced Chi2",800,800);
	c1->SetRightMargin(0.15);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);
	c1->SetTopMargin(0.08);

	hGridSearchReducedChi2->Draw("COLZ");
	hGridSearchReducedChi2->GetXaxis()->SetTitle("#phi_{low} (#circ)");
	hGridSearchReducedChi2->GetYaxis()->SetTitle("#phi_{high} (#circ)");
	hGridSearchReducedChi2->SetStats(0);

	//draw labels on each bin
	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){

			double val = hGridSearchReducedChi2->GetBinContent(ix,iy);
			if(val <= 0) continue;

			int globalbin = hGridSearchReducedChi2->GetBin(ix,iy);
			int rank = rankmap[globalbin];

			double x = hGridSearchReducedChi2->GetXaxis()->GetBinCenter(ix);
			double y = hGridSearchReducedChi2->GetYaxis()->GetBinCenter(iy);

			double dy = 0.1*(hGridSearchReducedChi2->GetYaxis()->GetBinWidth(iy));

			TString text = Form("#%d",rank);
			latex.DrawLatex(x,y+dy,text);

			text = Form("%.2f",val);
			latex.DrawLatex(x,y-dy,text);

		}
	}

	//force equal scaling on both x and y axis w/ gpad:
	gPad->SetFixedAspectRatio();
	outfile->cd();
	c1->Write("cGridSearchReducedChi2");
	c1->SaveAs("Lithium6_1plus_fourpixelchisquared.png");

	outfile->Close();

}

void Lithium6_1plus_fourpixelchisquared_2(){

	TString anglestring = "16782278";

	std::vector<TString> beamstringsx = {
											"gaus0025",
											"gaus003",
											"gaus0035",
											"gaus004",
											"gaus0045",
											"gaus005"
										};

	std::vector<TString> beamstringsy = {
											"gaus0005",
											"gaus001",
											"gaus0015",
											"gaus002",
											"gaus0025",
											"gaus003"
										};										

	double xEdges[7] = {0.00225, 0.00275, 0.00325, 0.00375, 0.00425, 0.00475, 0.00525};
	double yEdges[7] = {0.00025, 0.00075, 0.00125, 0.00175, 0.00225, 0.00275, 0.00325};										

	TH2D *hGridSearchChi2_2 = new TH2D("hGridSearchChi2_2", "GridSearchChi2_2", 6, xEdges, 6, yEdges);
	TH2D *hGridSearchReducedChi2_2 = new TH2D("hGridSearchReducedChi2_2", "GridSearchReducedChi2_2", 6, xEdges, 6, yEdges);

	//---------------------------------------------------
	//				Establish data values
	//---------------------------------------------------

	//uncomment for DESKTOP:
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
	//uncomment for LAPTOP:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";


	// TString path_pix_r71_w29 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r71_w29";
	// TString path_pix_r72_w29 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r72_w29";
	// TString path_pix_r72_w30 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r72_w30";
	// TString path_pix_r71_w30 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r71_w30";

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	// TH1 *hpix_r71_w29 = dynamic_cast<TH1*>(datafile->Get(path_pix_r71_w29));
	// TH1 *hpix_r72_w29 = dynamic_cast<TH1*>(datafile->Get(path_pix_r72_w29));
	// TH1 *hpix_r72_w30 = dynamic_cast<TH1*>(datafile->Get(path_pix_r72_w30));
	// TH1 *hpix_r71_w30 = dynamic_cast<TH1*>(datafile->Get(path_pix_r71_w30));
	// if(!hpix_r71_w29 || !hpix_r72_w29 || !hpix_r72_w30 || !hpix_r71_w30){
	// 	std::cerr << "Error retrieving at least one data histogram!" << std::endl;
	// 	return;
	// }

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram!" << std::endl;
		return;
	}

	// hpix_r71_w29->SetDirectory(0);
	// hpix_r72_w29->SetDirectory(0);
	// hpix_r72_w30->SetDirectory(0);
	// hpix_r71_w30->SetDirectory(0);
	hData->SetDirectory(0);
	datafile->Close();

	// double counts_pix_r71_w29 = hpix_r71_w29->GetEntries();
	// double counts_pix_r72_w29 = hpix_r72_w29->GetEntries();
	// double counts_pix_r72_w30 = hpix_r72_w30->GetEntries();
	// double counts_pix_r71_w30 = hpix_r71_w30->GetEntries();

	double dataIntegral = (hData->GetBinContent(hData->GetBin(30,8))) + (hData->GetBinContent(hData->GetBin(30,9))) + (hData->GetBinContent(hData->GetBin(31,9))) + (hData->GetBinContent(hData->GetBin(31,8)));

	// std::vector<double> fourpix_relcounts = {counts_pix_r71_w29/fourpixsum, counts_pix_r72_w29/fourpixsum, counts_pix_r72_w30/fourpixsum, counts_pix_r71_w30/fourpixsum};

	//for(const auto& fn : filenames){
	for(const auto& bsx : beamstringsx){
		for(const auto& bsy : beamstringsy){

			TString fn = Form("kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%sx_%sy_histos.root", anglestring.Data(), bsx.Data(), bsy.Data());

			//establish sim values:
			//uncomment for DESKTOP
			// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/%s",fn.Data());
			// TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
			//uncomment for LAPTOP:
			TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/%s",fn.Data());
			TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";


			TFile *simfile = new TFile(simFilePath,"READ");
			if(!simfile || simfile->IsZombie()){
				std::cerr << "Error opening sim file " << simFilePath << "\n\n";
				continue;
			}

			// TH1 *hpix_r71_w29_sim = dynamic_cast<TH1*>(simfile->GetPath(path_pix_r71_w29));
			// TH1 *hpix_r72_w29_sim = dynamic_cast<TH1*>(simfile->GetPath(path_pix_r72_w29));
			// TH1 *hpix_r72_w30_sim = dynamic_cast<TH1*>(simfile->GetPath(path_pix_r72_w30));
			// TH1 *hpix_r71_w30_sim = dynamic_cast<TH1*>(simfile->GetPath(path_pix_r71_w30));
			// if(!hpix_r71_w29_sim || !hpix_r72_w29_sim || !hpix_r72_w30_sim || !hpix_r71_w30_sim){
			// 	std::cerr << "Error retrieving at least one sim histogram!" << std::endl;
			// 	continue;
			// }

			TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
			if(!hSim){
				std::cerr << "Error retrieving sim histo from " << simFilePath << "\n\n";
				continue;
			}

			double simIntegral = (hSim->GetBinContent(hSim->GetBin(30,8))) + (hSim->GetBinContent(hSim->GetBin(30,9))) + (hSim->GetBinContent(hSim->GetBin(31,9))) + (hSim->GetBinContent(hSim->GetBin(31,8)));

			// hpix_r71_w29_sim->SetDirectory(0);
			// hpix_r72_w29_sim->SetDirectory(0);
			// hpix_r72_w30_sim->SetDirectory(0);
			// hpix_r71_w30_sim->SetDirectory(0);
			hSim->SetDirectory(0);
			simfile->Close();

			if( hSim->GetNbinsX() != hData->GetNbinsX() ){
				std::cerr << "Histogram binning does not match for data/sim histogram pair\n";
				continue;
			}

			// double simIntegral_r71_w29 = hpix_r71_w29_sim->Integral();
			// double simIntegral_r72_w29 = hpix_r72_w29_sim->Integral();
			// double simIntegral_r72_w30 = hpix_r72_w30_sim->Integral();
			// double simIntegral_r71_w30 = hpix_r71_w30_sim->Integral();

			// double scaleFactor_r71_w29;// = counts_pix_r71_w29/simIntegral_r71_w29;
			// double scaleFactor_r72_w29;// = counts_pix_r72_w29/simIntegral_r72_w29;
			// double scaleFactor_r72_w30;// = counts_pix_r72_w30/simIntegral_r72_w30;
			// double scaleFactor_r71_w30;// = counts_pix_r71_w30/simIntegral_r71_w30;

			// if(simIntegral_r71_w29!=0 && simIntegral_r72_w29!=0 && simIntegral_r72_w30!=0 && simIntegral_r71_w30!=0){

			// 	scaleFactor_r71_w29 = counts_pix_r71_w29/simIntegral_r71_w29;
			// 	scaleFactor_r72_w29 = counts_pix_r72_w29/simIntegral_r72_w29;
			// 	scaleFactor_r72_w30 = counts_pix_r72_w30/simIntegral_r72_w30;
			// 	scaleFactor_r71_w30 = counts_pix_r71_w30/simIntegral_r71_w30;

			// 	hpix_r71_w29_sim->Scale(scaleFactor_r71_w29);
			// 	hpix_r72_w29_sim->Scale(scaleFactor_r72_w29);
			// 	hpix_r72_w30_sim->Scale(scaleFactor_r72_w30);
			// 	hpix_r71_w30_sim->Scale(scaleFactor_r71_w30);

			// } else {
			// 	std::cerr << "sim histogram has zero integral, cannot be scaled...continuing...\n";
			// 	continue;
			// }



			double scaleFactor;
			if(simIntegral > 0){
				scaleFactor = dataIntegral/simIntegral;
				hSim->Scale(scaleFactor);
			} else {
				std::cerr << "sim histogram from " << fn.Data() << " has zero integral and thus cannot be scaled! Continuing..." << std::endl;
				continue;
			}

			// double counts_sim_pix_r71_w29 = hpix_r71_w29_sim->Integral();
			// double counts_sim_pix_r72_w29 = hpix_r72_w29_sim->Integral();
			// double counts_sim_pix_r72_w30 = hpix_r72_w30_sim->Integral();
			// double counts_sim_pix_r71_w30 = hpix_r71_w30_sim->Integral();

			// double fourpixsum_sim = counts_sim_pix_r71_w29 + counts_sim_pix_r72_w29 + counts_sim_pix_r72_w30 + counts_sim_pix_r71_w30;

			//std::vector<double> fourpix_relcounts_sim = {counts_sim_pix_r71_w29/fourpixsum_sim, counts_sim_pix_r72_w29/fourpixsum_sim, counts_sim_pix_r72_w30/fourpixsum_sim, counts_sim_pix_r71_w30/fourpixsum_sim};

			//calculate the chi squared for this sim file:

			double r71w29_data = (hData->GetBinContent(hData->GetBin(30,8)));
			double r72w29_data = (hData->GetBinContent(hData->GetBin(30,9)));
			double r72w30_data = (hData->GetBinContent(hData->GetBin(31,9)));
			double r71w30_data = (hData->GetBinContent(hData->GetBin(31,8)));

			double r71w29_data_error = (hData->GetBinError(hData->GetBin(30,8)));
			double r72w29_data_error = (hData->GetBinError(hData->GetBin(30,9)));
			double r72w30_data_error = (hData->GetBinError(hData->GetBin(31,9)));
			double r71w30_data_error = (hData->GetBinError(hData->GetBin(31,8)));

			double r71w29_sim = (hSim->GetBinContent(hSim->GetBin(30,8)));
			double r72w29_sim = (hSim->GetBinContent(hSim->GetBin(30,9)));
			double r72w30_sim = (hSim->GetBinContent(hSim->GetBin(31,9)));
			double r71w30_sim = (hSim->GetBinContent(hSim->GetBin(31,8)));

			double r71w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,8)));
			double r72w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,9)));
			double r72w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,9)));
			double r71w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,8)));

			double chi2 = ( std::pow((r71w29_data - r71w29_sim),2) / ( std::pow(r71w29_data_error,2) + std::pow(r71w29_sim_error,2)) ) + 
						  ( std::pow((r72w29_data - r72w29_sim),2) / ( std::pow(r72w29_data_error,2) + std::pow(r72w29_sim_error,2)) ) +
						  ( std::pow((r72w30_data - r72w30_sim),2) / ( std::pow(r72w30_data_error,2) + std::pow(r72w30_sim_error,2)) ) +
						  ( std::pow((r71w30_data - r71w30_sim),2) / ( std::pow(r71w30_data_error,2) + std::pow(r71w30_sim_error,2)) );

			Int_t ndf = 4 - 1;//4 pixels minus the 1 DOF (the scale factor is only fitted parameter here)

			double reducedchi2 = chi2/ndf;

			//fill chi2 histogram from grid search:

			std::map<TString,double> beamXmap = {{"gaus0025",0.0025},{"gaus003",0.003},{"gaus0035",0.0035},{"gaus004",0.004},{"gaus0045",0.0045},{"gaus005",0.005}};

			std::map<TString,double> beamYmap = {{"gaus0005",0.0005},{"gaus001",0.001},{"gaus0015",0.0015},{"gaus002",0.002},{"gaus0025",0.0025},{"gaus003",0.003}};

			double xbin = beamXmap[bsx];
			double ybin = beamYmap[bsy];

			hGridSearchChi2_2->Fill(xbin,ybin,chi2);
			hGridSearchReducedChi2_2->Fill(xbin,ybin,reducedchi2);
		}
	}

	TFile *outfile = new TFile("GridSearch2.root","RECREATE");

	outfile->cd();
	//hData->Write();
	//hSim->Write();
	hGridSearchChi2_2->Write();
	hGridSearchReducedChi2_2->Write();

	//set up square TCanvas here:
	TCanvas *c1 = new TCanvas("c1","Grid Search Reduced Chi2",800,800);
	c1->SetRightMargin(0.15);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);
	c1->SetTopMargin(0.08);

	hGridSearchReducedChi2_2->Draw("COLZ");
	hGridSearchReducedChi2_2->GetXaxis()->SetTitle("Beam Sigma X (m)");
	hGridSearchReducedChi2_2->GetYaxis()->SetTitle("Beam Sigma Y (m)");
	hGridSearchReducedChi2_2->SetStats(0);

	//force equal scaling on both x and y axis w/ gpad:
	gPad->SetFixedAspectRatio();
	outfile->cd();
	c1->Write("cGridSearchReducedChi2_2");
	c1->SaveAs("Lithium6_1plus_fourpixelchisquared_2.png");

	outfile->Close();//needs to be converted for tighter beamspot varying

}

void Lithium6_1plus_sixteenpixelchisquared(){

	TString anglestring = "178218";

	std::vector<TString> beamstrings = {
											//"fixed"
											"gaus001"
											// "gaus002",
											// "gaus003",
											// "gaus004",
											// "gaus005"
										};

	//std::vector<double> phis = {1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0};
	std::vector<double> phis = {2.0, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3.0, 3.125, 3.25, 3.375, 3.5};								
	std::vector<TString> phistrings;
	for(const auto& philower : phis){
		for(const auto& phiupper : phis){
			//phistrings.push_back(Form("%f_%f",-philower,phiupper));
			TString s = TString::Format("%s_%s", smartFormat(-philower).c_str(), smartFormat(phiupper).c_str());
			phistrings.push_back(s);
		}
	}

	// double xedges[8] = {-3.125, -2.875, -2.625, -2.375, -2.125, -1.875, -1.625, -1.375};
	// double yedges[8] = {1.375, 1.625, 1.875, 2.125, 2.375, 2.625, 2.875, 3.125};

	double xedges[14] = {-3.5625, -3.4375, -3.3125, -3.1875, -3.0625, -2.9375, -2.8125, -2.6875, -2.5625, -2.4375, -2.3125, -2.1875, -2.0625, -1.9375};
	double yedges[14] = {1.9375, 2.0625, 2.1875, 2.3125, 2.4375, 2.5625, 2.6875, 2.8125, 2.9375, 3.0625, 3.1875, 3.3125, 3.4375, 3.5625};

	TH2D *hGridSearchChi2 = new TH2D("hGridSearchChi2", "GridSearchChi2", 13, xedges, 13, yedges);
	TH2D *hGridSearchReducedChi2 = new TH2D("hGridSearchReducedChi2", "GridSearchReducedChi2", 13, xedges, 13, yedges);


	//---------------------------------------------------
	//				Establish data values
	//---------------------------------------------------

	//uncomment for DESKTOP:
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
	//uncomment for LAPTOP:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram!" << std::endl;
		return;
	}

	hData->SetDirectory(0);
	datafile->Close();

	//ring (y) range for 1+ groudn state from data
	//same with wedge/x range but this does not change between 1+ GS and 0+ 3.563 MeV state
	int x_minbin = 28+1;//1 added to account for root's first non-underflow bin is bin=1, which corresponds to ring/wedge 0 for us
	int x_maxbin = 31+1;
	int y_minbin = 6+1;
	int y_maxbin = 9+1;
	double dataIntegral = hData->Integral(x_minbin, x_maxbin, y_minbin, y_maxbin);

	for(const auto& ps : phistrings){
		for(const auto& bsx : beamstrings){
			for(const auto& bsy : beamstrings){

				TString fn = Form("kin2mc_7Li3He4He6Ligs_7500keV_theta%s_phi_%s_%sx_%sy_histos.root", anglestring.Data(), ps.Data(), bsx.Data(), bsy.Data());

				//uncomment for DESKTOP:
				//TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/%s",fn.Data());
				//TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

				//uncomment for LAPTOP:
				TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/%s",fn.Data());
				TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

				TFile *simfile = new TFile(simFilePath,"READ");
				if(!simfile || simfile->IsZombie()){
					std::cerr << "Error opening sim file " << simFilePath << "\n\n";
					continue;
				}

				TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
				if(!hSim){
					std::cerr << "Error retrieving sim histo from " << simFilePath << "\n\n";
					continue;
				}

				double simIntegral = hSim->Integral(x_minbin, x_maxbin, y_minbin, y_maxbin);

				hSim->SetDirectory(0);
				simfile->Close();

				if(hSim->GetNbinsX() != hData->GetNbinsX()){
					std::cerr << "Histogram binning does not match for data/sim histogram pair\n\n";
					continue;
				}

				double scaleFactor;
				if(simIntegral > 0){
					scaleFactor = dataIntegral/simIntegral;
					hSim->Scale(scaleFactor);
				} else {
					std::cerr << "sim histo from " << fn.Data() << " has zero integral and thus cannot be scaled!\n\n";
					continue;
				}

				double chi2 = 0.;
				double reducedchi2 = 0.;
				int ndf = 16 - 1;//16 pixels minus 1 DOF (scale factor only DOF on a per-run basis)

				for(int x=x_minbin; x<=x_maxbin; x++){
					
					for(int y=y_minbin; y<=y_maxbin; y++){

						double datavalue = hData->GetBinContent(hData->GetBin(x,y));
						double simvalue = hSim->GetBinContent(hSim->GetBin(x,y));

						if(datavalue == 0 && simvalue == 0){
							ndf -= 1;
							std::cout << "passing for x = " << x << " and y = " << y << "\n";
							continue;
						}
						
						double dataerror = hData->GetBinError(hData->GetBin(x,y));
						double simerror = hSim->GetBinError(hSim->GetBin(x,y));

						double numerator = std::pow((datavalue - simvalue) , 2);
						double denominator = std::pow(dataerror, 2) + std::pow(simerror, 2);
						chi2 += (numerator/denominator);

					}
				}

				reducedchi2 = chi2/ndf;

				int xbin = 0;
				int ybin = 0;

				// if(bsx == "fixed") xbin = 0;
				// else if(bsx == "gaus001") xbin = 1;
				// else if(bsx == "gaus002") xbin = 2;
				// else if(bsx == "gaus003") xbin = 3;
				// else if(bsx == "gaus004") xbin = 4;
				// else if(bsx == "gaus005") xbin = 5;

				// if(bsy == "fixed") ybin = 0;
				// else if(bsy == "gaus001") ybin = 1;
				// else if(bsy == "gaus002") ybin = 2;
				// else if(bsy == "gaus003") ybin = 3;
				// else if(bsy == "gaus004") ybin = 4;
				// else if(bsy == "gaus005") ybin = 5;

				std::pair<TString,TString> philowhigh = parsephistring(ps);
				TString philow = philowhigh.first;
				TString phihigh = philowhigh.second;

				hGridSearchChi2->Fill(philow.Atof(),phihigh.Atof(),chi2);
				hGridSearchReducedChi2->Fill(philow.Atof(),phihigh.Atof(),reducedchi2);
			}
		}
	}

	//add labels showing values and relative ranks on each bin:
	int nx = hGridSearchReducedChi2->GetNbinsX();
	int ny = hGridSearchReducedChi2->GetNbinsY();

	//store (value, global_bin) for ranking
	std::vector<std::pair<double,int>> values;
	values.reserve(nx*ny);

	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){
			double val = hGridSearchReducedChi2->GetBinContent(ix,iy);
			int globalbin = hGridSearchReducedChi2->GetBin(ix,iy);

			if(val > 0) values.push_back({val,globalbin});
		}
	}

	//sort by reduced chi2 value:
	std::sort(values.begin(), values.end(),
		      [](auto&a, auto&b){ return a.first < b.first; });

	std::map<int,int> rankmap;
	for(size_t i=0; i<values.size(); i++){
		rankmap[values[i].second] = i+1;
	}

	TLatex latex;
	latex.SetTextAlign(22);
	latex.SetTextSize(0.012);
	latex.SetTextColor(kBlack);

	TFile *outfile = new TFile("GridSearch_16pix.root","RECREATE");
	outfile->cd();
	hGridSearchChi2->Write();
	hGridSearchReducedChi2->Write();

	TCanvas *c1 = new TCanvas("c1","Grid Seach Reduced Chi2",800,800);
	c1->SetRightMargin(0.15);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);
	c1->SetTopMargin(0.08);

	hGridSearchReducedChi2->Draw("COLZ");
	hGridSearchReducedChi2->GetXaxis()->SetTitle("#phi_{low} (#circ)");
	hGridSearchReducedChi2->GetYaxis()->SetTitle("#phi_{high} (#circ)");
	hGridSearchReducedChi2->SetStats(0);

	//draw labels on each bin:
	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){

			double val = hGridSearchReducedChi2->GetBinContent(ix,iy);
			if(val <= 0) continue;

			int globalbin = hGridSearchReducedChi2->GetBin(ix,iy);
			int rank = rankmap[globalbin];

			double x = hGridSearchReducedChi2->GetXaxis()->GetBinCenter(ix);
			double y = hGridSearchReducedChi2->GetYaxis()->GetBinCenter(iy);

			double dy = 0.1*(hGridSearchReducedChi2->GetYaxis()->GetBinWidth(iy));

			TString text = Form("#%d",rank);
			latex.DrawLatex(x,y+dy,text);

			text = Form("%.2f",val);
			latex.DrawLatex(x,y-dy,text);

		}
	}

	gPad->SetFixedAspectRatio();
	outfile->cd();
	c1->Write("cGridSearchReducedChi2");
	c1->SaveAs("Lithium6_1plus_sixteenpixelchisquared.png");

	outfile->Close();
}

void Lithium6_1plus_sixteenpixelchisquared_2(){

	TString anglestring = "16782278";

	std::vector<TString> beamstringsx = {
											"fixed",
											"gaus0005",
											"gaus001",
											"gaus0015",
											"gaus002",
											"gaus0025"
										};

	std::vector<TString> beamstringsy = {
											"fixed",
											"gaus0005",
											"gaus001",
											"gaus0015"										
										};

	//double xEdges[7] = {0.00225, 0.00275, 0.00325, 0.00375, 0.00425, 0.00475, 0.00525};
	//double yEdges[7] = {0.00025, 0.00075, 0.00125, 0.00175, 0.00225, 0.00275, 0.00325};

	// double xEdges[6] = {-0.00025, 0.00025, 0.00075, 0.00125, 0.00175, 0.00225};
	// double yEdges[6] = {-0.00025, 0.00025, 0.00075, 0.00125, 0.00175, 0.00225};

	// double xEdges[5] = {-0.000125, 0.000125, 0.000375, 0.000625, 0.000875};
	// double yEdges[5] = {-0.000125, 0.000125, 0.000375, 0.000625, 0.000875};

	double xEdges[7] = {-0.00025, 0.00025, 0.00075, 0.00125, 0.00175, 0.00225, 0.00275};
	double yEdges[5] = {-0.00025, 0.00025, 0.00075, 0.00125, 0.00175};

	TH2D *hGridSearchChi2_2 = new TH2D("hGridSearchChi2_2", "GridSearchChi2_2", 6, xEdges, 4, yEdges);
	TH2D *hGridSearchReducedChi2_2 = new TH2D("hGridSearchReducedChi2_2", "GridSearchReducedChi2_2", 6, xEdges, 4, yEdges);

	//---------------------------------------------------
	//				Establish data values
	//---------------------------------------------------

	//uncomment for DESKTOP:
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
	//uncomment for LAPTOP:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram!" << std::endl;
		return;
	}

	hData->SetDirectory(0);
	datafile->Close();

	int x_minbin = 29;
	int x_maxbin = 32;
	int y_minbin = 7;
	int y_maxbin = 10;
	double dataIntegral = hData->Integral(x_minbin, x_maxbin, y_minbin, y_maxbin);

	for(const auto& bsx : beamstringsx){
		for(const auto& bsy : beamstringsy){

			TString fn = Form("kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%sx_%sy_histos.root",anglestring.Data(),bsx.Data(),bsy.Data());

			//uncomment for DESKTOP:
			//TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/%s",fn.Data());
			//TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

			//uncomment for LAPTOP:
			TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/%s",fn.Data());
			TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

			TFile *simfile = new TFile(simFilePath,"READ");
			if(!simfile || simfile->IsZombie()){
				std::cerr << "Error opening sim file " << simFilePath << "\n\n";
				continue;
			}

			TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
			if(!hSim){
				std::cerr << "Error retrieving sim histo from " << simFilePath << "\n\n";
				continue;
			}

			double simIntegral = hSim->Integral(x_minbin, x_maxbin, y_minbin, y_maxbin);

			hSim->SetDirectory(0);
			simfile->Close();

			if(hSim->GetNbinsX() != hData->GetNbinsX()){
				std::cerr << "Histogram binning does not match for data/sim histogram pair\n\n";
				continue;
			}

			double scaleFactor;
			if(simIntegral > 0){
				scaleFactor = dataIntegral/simIntegral;
				hSim->Scale(scaleFactor);
			} else {
				std::cerr << "sim histo from " << fn.Data() << " has zero integral and thus cannot be scaled!\n\n";
				continue;
			}

			double chi2 = 0.;
			double reducedchi2 = 0.;
			int ndf = -666;

			for(int x=x_minbin; x<=x_maxbin; x++){
				ndf = 16 - 1;
				for(int y=y_minbin; y<=y_maxbin; y++){

					double datavalue = hData->GetBinContent(hData->GetBin(x,y));
					double simvalue = hSim->GetBinContent(hSim->GetBin(x,y));

					if(datavalue == 0 && simvalue == 0){
						ndf -= 1;
						std::cout << "passing for x = " << x << " and y = " << y << "\n";
						continue;
					}

					double dataerror = hData->GetBinError(hData->GetBin(x,y));
					double simerror = hSim->GetBinError(hSim->GetBin(x,y));

					double numerator = std::pow((datavalue - simvalue), 2);
					double denominator = std::pow(dataerror, 2) + std::pow(simerror, 2);
					chi2 += (numerator/denominator);

				}
			}

			reducedchi2 = chi2/ndf;

			double xbin = -1;
			double ybin = -1;

			if(bsx == "fixed") xbin = 0;
			else if(bsx == "gaus0005") xbin = 0.0005;
			else if(bsx == "gaus001") xbin = 0.001;
			else if(bsx == "gaus0015") xbin = 0.0015;
			else if(bsx == "gaus002") xbin = 0.002;
			else if(bsx == "gaus0025") xbin = 0.0025;

			if(bsy == "fixed") ybin = 0;
			else if(bsy == "gaus0005") ybin = 0.0005;
			else if(bsy == "gaus001") ybin = 0.001;
			else if(bsy == "gaus0015") ybin = 0.0015;

			hGridSearchChi2_2->Fill(xbin, ybin, chi2);
			hGridSearchReducedChi2_2->Fill(xbin, ybin, reducedchi2);
		}
	}

	TFile *outfile = new TFile("GridSearch_16pix_2.root","RECREATE");
	outfile->cd();
	hGridSearchChi2_2->Write();
	hGridSearchReducedChi2_2->Write();

	TCanvas *c1 = new TCanvas("c1","Grid Seach Reduced Chi2",800,800);
	c1->SetRightMargin(0.15);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);
	c1->SetTopMargin(0.08);

	hGridSearchReducedChi2_2->Draw("COLZ");
	hGridSearchReducedChi2_2->GetXaxis()->SetTitle("Beam Sigma X (m)");
	hGridSearchReducedChi2_2->GetYaxis()->SetTitle("Beam Sigma Y (m)");
	hGridSearchReducedChi2_2->SetStats(0);

	gPad->SetFixedAspectRatio();
	outfile->cd();
	c1->Write("cGridSearchReducedChi2_2");
	c1->SaveAs("Lithium6_1plus_sixteenpixelchisquared_2.png");

	outfile->Close();
}



void Lithium6_1plus_chiSquaredFromFourPixelsEnergySpectra(){
	//this function calculates a cumulative chi^2 for comparing simulation data with experimental data
	//The simulation data is scaled to match experimental data by using total events
	//This scaled simulation data is then passed to the experimental data for the ROOT TH1::Chi2Test()
	//This is done for the four main pixels illuminated by the reaction: 
	//								(71,29), (72,29), (72,30), (71,30)

	const int rebinfactor = 4;

	std::vector<std::pair<int,int>> ringwedgepairs = {{71,29}, {72,29}, {72,30}, {71,30}};

	// std::vector<TString> anglestrings = {
	// 										"188208",
	// 										"178218",
	// 										"168228",
	// 										"158238",
	// 										"148248"
	// 									};

	std::vector<std::pair<TString, int>> anglestrings_binnums = {
																	{"188208",1},
																	{"178218",2},
																	{"168228",3},
																	{"158238",4},
																	{"148248",5}
																 };

	// std::vector<TString> beamstringsx = {
	// 										"fixed"
	// 									};

	// std::vector<TString> beamstringsy = {
	// 										"fixed"
	// 									};

	TString beamstringx = "gaus001";
	TString beamstringy = "gaus001";

	int SABRE_ID = 3;

	TFile *outfile = new TFile("EjectileAngleScan.root","RECREATE");		


	TH1D *hChi2 = new TH1D("hChi2","Chi2", 5, 0.5, 5.5);
	TH1D *hReducedChi2 = new TH1D("hReducedChi2","ReducedChi2", 5, 0.5, 5.5);

	//uncomment below for DESKTOP:
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = Form("SABRE/SABRE%d/Pixels/",SABRE_ID);
	// dataHistLocalPath = dataHistLocalPath + pixelhistoname;

	//---------------------------------------------------
	//				Establish data file
	//---------------------------------------------------

	//uncomment below for LAPTOP:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		std::cerr << "dataFilePath = " << dataFilePath << "\n\n";
		return;// "ERROR";
	}


	//loop through the angle strings, as each angle string is an entry in the chi2 histogram:
	for(const auto& angstr_bn : anglestrings_binnums){

		//establish data holders for chi2, p, and ndf of each pixel:
		std::vector<std::array<double,3>> chi2_results;

		for(const auto& rwpair : ringwedgepairs){

			int ring = rwpair.first;
			int wedge = rwpair.second;

			TString pixelhistoname = Form("SABRE/SABRE%d/Pixels/hSABRE%d_pixel_r%dw%d_ESummary",SABRE_ID,SABRE_ID,ring,wedge);

			//---------------------------------------------------
			//				Retrieve Data Histogram
			//---------------------------------------------------
				TH1 *hData = dynamic_cast<TH1*>(datafile->Get(pixelhistoname));
				if(!hData){
					std::cerr << "Error retrieving data histogram" << std::endl;
					std::cerr << "dataFilePath = " << dataFilePath << "\n";
					//std::cerr << "dataHistLocalPath = " << dataHistLocalPath << "\n\n";
					datafile->Close();
					return;// "ERROR";
				}
				hData->SetDirectory(0);

			//---------------------------------------------------
			//				Establish sim values
			//---------------------------------------------------

			//uncomment below for DESKTOP:
			// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%sx_%sy_histos.root",angstr_bn.first.Data(),beamstringx.Data(),beamstringy.Data());
			// TString simHistLocalPath = Form("SABRE/SABRE%d/Pixels/",SABRE_ID);
			// simHistLocalPath = simHistLocalPath + pixelhistoname;

			//uncomment below for LAPTOP:
			TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%sx_%sy_histos.root",angstr_bn.first.Data(),beamstringx.Data(),beamstringy.Data());

			TFile *simfile = new TFile(simFilePath,"READ");
			if(!simfile || simfile->IsZombie()){
				std::cerr << "Error opening sim file" << std::endl;
				return;// "ERROR";
			}

			TH1 *hSim = dynamic_cast<TH1*>(simfile->Get(pixelhistoname));
			if(!hSim){
				std::cerr << "Error retrieving sim histogram" << std::endl;
				simfile->Close();
				return;// "ERROR";
			}
			hSim->SetDirectory(0);
			simfile->Close();


			//---------------------------------------------------
			//					 Check Binning
			//---------------------------------------------------

			if((hData->GetNbinsX() != hSim->GetNbinsX()) || (hData->GetNbinsY() != hSim->GetNbinsY())){
				std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
				return;// "ERROR";
			}

			//---------------------------------------------------
			//				Determine Scale Factor
			//					  (and scale)
			//---------------------------------------------------

				double dataIntegral = hData->Integral();
				double simIntegral = hSim->Integral();

				double scaleFactor = 0.;

				if(simIntegral > 0){
					scaleFactor = dataIntegral/simIntegral;
					//hSim->Scale(scaleFactor); //scaling hsim in chi2 function now
					//std::cout << "Applied global scale factor to sim: " << scaleFactor << std::endl;
				} else {
					std::cerr << "sim histogram has zero integral and thus cannot scale" << std::endl;
					return;// "ERROR";
				}

			//---------------------------------------------------
			//			   Rebin Sim and Data Histos
			//---------------------------------------------------
				//rebinfactor defined at top of function, typically = 4

				//std::cout << "Rebinning..." << std::endl;
				hSim->Rebin(rebinfactor);
				hData->Rebin(rebinfactor);
				// hSim->Delete();
				// hData->Delete();
				//std::cout << "Done!" << std::endl;


			//---------------------------------------------------
			//			  Calculate and Store Chi^2
			//---------------------------------------------------

				// std::cout << "calculating chi2..." << std::endl;
				// double results[3];
				// //std::array<double,3> results;
				// std::cout << "results declared..." << std::endl;
				// //hData_rebin->Chi2Test(hSim_rebin, "CHI2 UW", results);
				// double chi2 = hData->Chi2Test(hSim, "UW CHI2");
				// double ndf = hData->Chi2Test(hSim, "UW NDF");
				// double pval = hData->Chi2Test(hSim, "UW p");
				// results[0] = chi2;
				// results[1] = ndf;
				// results[2] = pval;
				// std::cout << "results now holds chi2 results!" << std::endl;
				// //results[0] = chi2
				// //results[1] = ndf
				// //results[2] = pval

				// chi2_results.push_back({ results[0], results[1], results[2] });
				// std::cout << "done!" << std::endl;

				//we can't trust whatever is going on int Chi2Test, so let's just write our own helper function elsewhere
				double ndf;
				//double chi2 = ComputeChi2BinByBin(hData, hSim, ndf);
				double chi2 = ComputeChi2BinByBin(hData, hSim, scaleFactor, ndf);
				double pval = TMath::Prob(chi2,ndf);

				std::array<double,3> results = {chi2, ndf, pval};
				chi2_results.push_back(results);

			//---------------------------------------------------
			//			  		Delete Histograms
			//---------------------------------------------------

				std::cout << "deleting..." << std::endl;
				outfile->cd();
				// hSim->Write();
				// hData->Write();
				delete hSim;
				delete hData;
				// delete hSim_rebin;
				// delete hData_rebin;
				std::cout << "deleted!" << std::endl;

		}

		//at this point, we have all (ring,wedge) pairs done for a given angstr
		//and the results for each pixel are appended in the chi2_results vector!
		//recall that:
		//	results[0] = chi2
		//	results[1] = ndf
		//	results[2] = pval

		//we know that the cumulative chi^2 is just the sum of all the chi^2:
		double totalchi2 = 0;
		double totalndf = 0;
		for(const auto& res : chi2_results){
			totalchi2 += res[0];
			totalndf += res[1];
		}
		double totalreducedchi2 = totalchi2/totalndf;

		//fill histograms for this angstr bin:
		int bin = angstr_bn.second;

		// hChi2->Fill(bin, totalchi2);
		// hReducedChi2->Fill(bin, totalreducedchi2);

		hChi2->SetBinContent(bin, totalchi2);
		hReducedChi2->SetBinContent(bin, totalreducedchi2);

	}

	datafile->Close();

	std::vector<TString> labels = {"19.8#circ #pm 1#circ",
								   "19.8#circ #pm 2#circ",
								   "19.8#circ #pm 3#circ",
								   "19.8#circ #pm 4#circ",
								   "19.8#circ #pm 5#circ"};

	for(size_t i = 0; i < labels.size(); i++){
		hChi2->GetXaxis()->SetBinLabel(i+1, labels[i]);
		//hChi2->LabelsOption("v");

		hReducedChi2->GetXaxis()->SetBinLabel(i+1, labels[i]);
		hReducedChi2->GetYaxis()->SetRangeUser(0,100);
		//hReducedChi2->LabelsOption("v");
	}


	//now we just need to save the histograms to a root file:
	outfile->cd();
	hChi2->Write();
	hReducedChi2->Write();


	//also just draw the reduced chi2 to screen for quick testing:
	// hReducedChi2->SetDirectory(0);
	// hReducedChi2->Draw();


	outfile->Close();
}

void Lithium6_1plus_BeamSpotSearch(){

	/*
		for a single theta range and a single phi range, search through beamspot profiles in a grid
		and produce histograms that show the chi^2 and the reduced chi^2 for the following criteria:

			i)		4pixel chi^2 (and chi^2/NDF)
			ii)		4+12(16)pixel chi^2 (and chi^2/NDF)
			iii)	4pixel-by-pixel energy spectrum chi^2 (and chi^2/NDF)
		

	*/

	const int rebinfactor = 4;//used to rebin histograms to more coarse spacing, helps for comparison

	int SABRE_ID = 3;

	std::vector<std::pair<int,int>> ringwedgepairs = {{71,29}, {72,29}, {72,30}, {71,30}};
	
	TString anglestring = "178218";
	TString phistring = "-2.125_2.125";
	std::vector<TString> beamstrings = {"gaus0005", "gaus001", "gaus0015", "gaus002", "gaus0025"};

//							   0      0.0005   0.001   0.0015   0.002    0.0025
	double xEdges[7] = {-0.00025, 0.00025, 0.00075, 0.00125, 0.00175, 0.00225, 0.00275};
	double yEdges[7] = {-0.00025, 0.00025, 0.00075, 0.00125, 0.00175, 0.00225, 0.00275};

	//--------------------------------------establish necessary histograms--------------------------------------
	//4pixel count histograms:
	TH2D *hBeamSpotSearch_4pixel_chi2 = new TH2D("hBeamSpotSearch_4pixel_chi2","4pixel #chi^{2}", 6, xEdges, 6, yEdges);
	TH2D *hBeamSpotSearch_4pixel_redchi2 = new TH2D("hBeamSpotSearch_4pixel_redchi2","BeamSpotSearch_4pixel_redchi2", 6, xEdges, 6, yEdges);

	//4+12(16)pixel count histograms:
	TH2D *hBeamSpotSearch_16pixel_chi2 = new TH2D("hBeamSpotSearch_16pixel_chi2", "BeamSpotSearch_16pixel_chi2", 6, xEdges, 6, yEdges);
	TH2D *hBeamSpotSearch_16pixel_redchi2 = new TH2D("hBeamSpotSearch_16pixel_redchi2", "BeamSpotSearch_16pixel_redchi2", 6, xEdges, 6, yEdges);

	//4pixel-by-pixel energy spectrum histograms
	TH2D *hBeamSpotSearch_4pixelspectra_chi2 = new TH2D("hBeamSpotSearch_4pixelspectra_chi2", "BeamSpotSearch_4pixelspectra_chi2", 6, xEdges, 6, yEdges);
	TH2D *hBeamSpotSearch_4pixelspectra_redchi2 = new TH2D("hBeamSpotSearch_4pixelspectra_redchi2", "BeamSpotSearch_4pixelspectra_redchi2", 6, xEdges, 6, yEdges);

	//-----------------------------------4, 4+12 pixel count histos ------------------------------------

	//------------------------------------retrieve data histograms--------------------------------------
	//first, let's open the data file and retrieve relevant histograms for 4pixel and 4+12(16)pixel calculations:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_1plus_output.root";
	TFile *datafile = new TFile(dataFilePath, "READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file: " << dataFilePath << std::endl;
		return;
	}

	TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
	TH2 *hData_RingsVSWedges = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData_RingsVSWedges){
		std::cerr << "Error retrieving data histogram " << dataHistLocalPath << " from file " << dataFilePath << std::endl;
		return;
	}
	hData_RingsVSWedges->SetDirectory(0);
	double dataIntegral = (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(30,8))) + (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(30,9))) + (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(31,9))) + (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(31,8)));

	int x_minbin = 28+1;//1 added to account for root's first non-underflow bin is bin=1, which corresponds to ring/wedge 0 for us
	int x_maxbin = 31+1;//1 added to account for root's first non-underflow bin is bin=1, which corresponds to ring/wedge 0 for us
	int y_minbin = 6+1;//1 added to account for root's first non-underflow bin is bin=1, which corresponds to ring/wedge 0 for us
	int y_maxbin = 9+1;//1 added to account for root's first non-underflow bin is bin=1, which corresponds to ring/wedge 0 for us
	double dataIntegral_16 = hData_RingsVSWedges->Integral(x_minbin, x_maxbin, y_minbin, y_maxbin);

	for(const auto& bsx : beamstrings){
		for(const auto& bsy : beamstrings){

			TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_phi_%s_%sx_%sy_histos.root", anglestring.Data(), phistring.Data(), bsx.Data(), bsy.Data());
			TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

			TFile *simfile = new TFile(simFilePath, "READ");
			if(!simfile || simfile->IsZombie()){
				std::cerr << "Error opening sim file " << simFilePath << "\n\n";
				continue;
			}

			TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
			if(!hSim){
				std::cerr << "Error retrieving sim histo from " << simFilePath << "\n\n";
				continue;
			}

			double simIntegral = (hSim->GetBinContent(hSim->GetBin(30,8))) + (hSim->GetBinContent(hSim->GetBin(30,9))) + (hSim->GetBinContent(hSim->GetBin(31,9))) + (hSim->GetBinContent(hSim->GetBin(31,8)));

			hSim->SetDirectory(0);
			

			if(hSim->GetNbinsX() != hData_RingsVSWedges->GetNbinsX()){
				std::cerr << "Histogram binning does not match for data/sim histogram pair\n";
				continue;
			}

			double scaleFactor;
			if(simIntegral > 0){
				scaleFactor = dataIntegral/simIntegral;
				hSim->Scale(scaleFactor);
			} else {
				std::cerr << "sim histogram from " << simFilePath.Data() << " has zero integral and thus cannot be scaled! Continuing..." << std::endl;
				continue;
			}

			//calculate the chi squared for this sim file:

			double r71w29_data = (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(30,8)));//binx = 30 is wedge 29, biny = 8 is local ring 7 which for SABRE3 is ring 71
			double r72w29_data = (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(30,9)));//binx = 30 is wedge 29, biny = 9 is local ring 8 which for SABRE3 is ring 72
			double r72w30_data = (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(31,9)));//binx = 31 is wedge 30, biny = 9 is local ring 8 which for SABRE3 is ring 72
			double r71w30_data = (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(31,8)));//binx = 31 is wedge 30, biny = 8 is local ring 7 which for SABRE3 is ring 71

			double r71w29_data_error = (hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(30,8)));
			double r72w29_data_error = (hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(30,9)));
			double r72w30_data_error = (hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(31,9)));
			double r71w30_data_error = (hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(31,8)));

			double r71w29_sim = (hSim->GetBinContent(hSim->GetBin(30,8)));
			double r72w29_sim = (hSim->GetBinContent(hSim->GetBin(30,9)));
			double r72w30_sim = (hSim->GetBinContent(hSim->GetBin(31,9)));
			double r71w30_sim = (hSim->GetBinContent(hSim->GetBin(31,8)));

			double r71w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,8)));
			double r72w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,9)));
			double r72w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,9)));
			double r71w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,8)));

			double chi2 = ( std::pow((r71w29_data - r71w29_sim),2) / ( std::pow(r71w29_data_error,2) + std::pow(r71w29_sim_error,2)) ) + 
				  		  ( std::pow((r72w29_data - r72w29_sim),2) / ( std::pow(r72w29_data_error,2) + std::pow(r72w29_sim_error,2)) ) +
						  ( std::pow((r72w30_data - r72w30_sim),2) / ( std::pow(r72w30_data_error,2) + std::pow(r72w30_sim_error,2)) ) +
						  ( std::pow((r71w30_data - r71w30_sim),2) / ( std::pow(r71w30_data_error,2) + std::pow(r71w30_sim_error,2)) );

			double ndf = 4 - 1;//4 pixels minus the 1 DOF (the scale factor is only fitted parameter here)

			double reducedchi2 = chi2/ndf;

			//fill chi2 histogram from grid search:
			double xbin = -1, ybin = -1;

			if(bsx == "gaus0005") xbin = 0.0005;
			else if(bsx == "gaus001") xbin = 0.001;
			else if(bsx == "gaus0015") xbin = 0.0015;
			else if(bsx == "gaus002") xbin = 0.002;
			else if(bsx == "gaus0025") xbin = 0.0025;

			if(bsy == "gaus0005") ybin = 0.0005;
			else if(bsy == "gaus001") ybin = 0.001;
			else if(bsy == "gaus0015") ybin = 0.0015;
			else if(bsy == "gaus002") ybin = 0.002;
			else if(bsy == "gaus0025") ybin = 0.0025;

			hBeamSpotSearch_4pixel_chi2->Fill(xbin, ybin, chi2);
			hBeamSpotSearch_4pixel_redchi2->Fill(xbin, ybin, reducedchi2);

			//now lets do the 4+12(16) pixels
			chi2 = 0.;
			reducedchi2 = 0.;
			ndf = 16-1;
			for(int x=x_minbin; x<=x_maxbin; x++){
				for(int y=y_minbin; y<=y_maxbin; y++){

					int datavalue = hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(x,y));
					int simvalue = hSim->GetBinContent(hSim->GetBin(x,y));

					if(datavalue == 0 && simvalue == 0){
						ndf -= 1;
						std::cout << "Passing for x = " << x << " and y = " << y << ". Decreasing NDF by 1 as result." << std::endl;
						continue;
					}

					double dataerror = hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(x,y));
					double simerror = hSim->GetBinError(hSim->GetBin(x,y));

					double numerator = std::pow((datavalue - simvalue), 2);
					double denominator = std::pow(dataerror, 2) + std::pow(simerror, 2);
					chi2 += (numerator/denominator);

				}
			}

			reducedchi2 = chi2/ndf;

			hBeamSpotSearch_16pixel_chi2->Fill(xbin, ybin, chi2);
			hBeamSpotSearch_16pixel_redchi2->Fill(xbin, ybin, reducedchi2);

			//now, let's get the pixel energy spectra necessary
			std::vector<std::array<double,3>> chi2_results;
			for(const auto& pixel : ringwedgepairs){

				int ring = pixel.first;
				int wedge = pixel.second;

				TString pixelhistoname = Form("SABRE/SABRE%d/Pixels/hSABRE%d_pixel_r%dw%d_ESummary", SABRE_ID, SABRE_ID, ring, wedge);

				TH1 *hData = dynamic_cast<TH1*>(datafile->Get(pixelhistoname));
				if(!hData){
					std::cerr << "Error retrieving data histogram" << std::endl;
					std::cerr << "dataFilePath = " << dataFilePath << "\n";
					return;
				}

				TH1 *hSim_ = dynamic_cast<TH1*>(simfile->Get(pixelhistoname));
				if(!hSim_){
					std::cerr << "Error retrieving sim histogram " << pixelhistoname << std::endl;
					std::cerr << "simFilePath = " << simFilePath << std::endl;
					return;
				}

				hSim_->SetDirectory(0);

				if((hData->GetNbinsX() != hSim_->GetNbinsX())){
					std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
					return;
				}

				dataIntegral = hData->Integral();
				simIntegral = hSim_->Integral();

				scaleFactor = 0.;

				if(simIntegral > 0){
					scaleFactor = dataIntegral/simIntegral;
				} else {
					std::cerr << "sim histogram has zero integral and thus cannot scale" << std::endl;
					return;
				}

				hSim_->Rebin(rebinfactor);
				hData->Rebin(rebinfactor);

				ndf = 0;
				chi2 = ComputeChi2BinByBin(hData, hSim_, scaleFactor, ndf);
				double pval = TMath::Prob(chi2,ndf);
				std::array<double,3> results = {chi2, ndf, pval};
				chi2_results.push_back(results);

				delete hSim_;
				hData->SetDirectory(0);
			}
			simfile->Close();

			double totalchi2 = 0.;
			double totalndf = 0.;
			for(const auto& res : chi2_results){
				totalchi2 += res[0];
				totalndf += res[1];
			}
			double totalreducedchi2 = totalchi2/totalndf;

			hBeamSpotSearch_4pixelspectra_chi2->Fill(xbin, ybin, totalchi2);
			hBeamSpotSearch_4pixelspectra_redchi2->Fill(xbin, ybin, totalreducedchi2);
		}
	}

	//do fixedx_fixedy here:
	beamstrings = {"fixed"};
	for(const auto& bsx : beamstrings){
		for(const auto& bsy : beamstrings){

			TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_phi_%s_%sx_%sy_histos.root", anglestring.Data(), phistring.Data(), bsx.Data(), bsy.Data());
			TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

			TFile *simfile = new TFile(simFilePath, "READ");
			if(!simfile || simfile->IsZombie()){
				std::cerr << "Error opening sim file " << simFilePath << "\n\n";
				continue;
			}

			TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
			if(!hSim){
				std::cerr << "Error retrieving sim histo from " << simFilePath << "\n\n";
				continue;
			}

			double simIntegral = (hSim->GetBinContent(hSim->GetBin(30,8))) + (hSim->GetBinContent(hSim->GetBin(30,9))) + (hSim->GetBinContent(hSim->GetBin(31,9))) + (hSim->GetBinContent(hSim->GetBin(31,8)));

			hSim->SetDirectory(0);
			

			if(hSim->GetNbinsX() != hData_RingsVSWedges->GetNbinsX()){
				std::cerr << "Histogram binning does not match for data/sim histogram pair\n";
				continue;
			}

			double scaleFactor;
			if(simIntegral > 0){
				scaleFactor = dataIntegral/simIntegral;
				hSim->Scale(scaleFactor);
			} else {
				std::cerr << "sim histogram from " << simFilePath.Data() << " has zero integral and thus cannot be scaled! Continuing..." << std::endl;
				continue;
			}

			//calculate the chi squared for this sim file:

			double r71w29_data = (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(30,8)));//binx = 30 is wedge 29, biny = 8 is local ring 7 which for SABRE3 is ring 71
			double r72w29_data = (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(30,9)));//binx = 30 is wedge 29, biny = 9 is local ring 8 which for SABRE3 is ring 72
			double r72w30_data = (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(31,9)));//binx = 31 is wedge 30, biny = 9 is local ring 8 which for SABRE3 is ring 72
			double r71w30_data = (hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(31,8)));//binx = 31 is wedge 30, biny = 8 is local ring 7 which for SABRE3 is ring 71

			double r71w29_data_error = (hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(30,8)));
			double r72w29_data_error = (hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(30,9)));
			double r72w30_data_error = (hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(31,9)));
			double r71w30_data_error = (hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(31,8)));

			double r71w29_sim = (hSim->GetBinContent(hSim->GetBin(30,8)));
			double r72w29_sim = (hSim->GetBinContent(hSim->GetBin(30,9)));
			double r72w30_sim = (hSim->GetBinContent(hSim->GetBin(31,9)));
			double r71w30_sim = (hSim->GetBinContent(hSim->GetBin(31,8)));

			double r71w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,8)));
			double r72w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,9)));
			double r72w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,9)));
			double r71w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,8)));

			double chi2 = ( std::pow((r71w29_data - r71w29_sim),2) / ( std::pow(r71w29_data_error,2) + std::pow(r71w29_sim_error,2)) ) + 
				  		  ( std::pow((r72w29_data - r72w29_sim),2) / ( std::pow(r72w29_data_error,2) + std::pow(r72w29_sim_error,2)) ) +
						  ( std::pow((r72w30_data - r72w30_sim),2) / ( std::pow(r72w30_data_error,2) + std::pow(r72w30_sim_error,2)) ) +
						  ( std::pow((r71w30_data - r71w30_sim),2) / ( std::pow(r71w30_data_error,2) + std::pow(r71w30_sim_error,2)) );

			double ndf = 4 - 1;//4 pixels minus the 1 DOF (the scale factor is only fitted parameter here)

			double reducedchi2 = chi2/ndf;

			//fill chi2 histogram from grid search:
			double xbin = 0, ybin = 0;

			// if(bsx == "gaus0005") xbin = 0.0005;
			// else if(bsx == "gaus001") xbin = 0.001;
			// else if(bsx == "gaus0015") xbin = 0.0015;
			// else if(bsx == "gaus002") xbin = 0.002;
			// else if(bsx == "gaus0025") xbin = 0.0025;

			// if(bsy == "gaus0005") ybin = 0.0005;
			// else if(bsy == "gaus001") ybin = 0.001;
			// else if(bsy == "gaus0015") ybin = 0.0015;
			// else if(bsy == "gaus002") ybin = 0.002;
			// else if(bsy == "gaus0025") ybin = 0.0025;

			hBeamSpotSearch_4pixel_chi2->Fill(xbin, ybin, chi2);
			hBeamSpotSearch_4pixel_redchi2->Fill(xbin, ybin, reducedchi2);

			//now lets do the 4+12(16) pixels
			chi2 = 0.;
			reducedchi2 = 0.;
			ndf = 16-1;
			for(int x=x_minbin; x<=x_maxbin; x++){
				for(int y=y_minbin; y<=y_maxbin; y++){

					int datavalue = hData_RingsVSWedges->GetBinContent(hData_RingsVSWedges->GetBin(x,y));
					int simvalue = hSim->GetBinContent(hSim->GetBin(x,y));

					if(datavalue == 0 && simvalue == 0){
						ndf -= 1;
						std::cout << "Passing for x = " << x << " and y = " << y << ". Decreasing NDF by 1 as result." << std::endl;
						continue;
					}

					double dataerror = hData_RingsVSWedges->GetBinError(hData_RingsVSWedges->GetBin(x,y));
					double simerror = hSim->GetBinError(hSim->GetBin(x,y));

					double numerator = std::pow((datavalue - simvalue), 2);
					double denominator = std::pow(dataerror, 2) + std::pow(simerror, 2);
					chi2 += (numerator/denominator);

				}
			}

			reducedchi2 = chi2/ndf;

			hBeamSpotSearch_16pixel_chi2->Fill(xbin, ybin, chi2);
			hBeamSpotSearch_16pixel_redchi2->Fill(xbin, ybin, reducedchi2);

			//now, let's get the pixel energy spectra necessary
			std::vector<std::array<double,3>> chi2_results;
			for(const auto& pixel : ringwedgepairs){

				int ring = pixel.first;
				int wedge = pixel.second;

				TString pixelhistoname = Form("SABRE/SABRE%d/Pixels/hSABRE%d_pixel_r%dw%d_ESummary", SABRE_ID, SABRE_ID, ring, wedge);

				TH1 *hData = dynamic_cast<TH1*>(datafile->Get(pixelhistoname));
				if(!hData){
					std::cerr << "Error retrieving data histogram" << std::endl;
					std::cerr << "dataFilePath = " << dataFilePath << "\n";
					return;
				}

				TH1 *hSim_ = dynamic_cast<TH1*>(simfile->Get(pixelhistoname));
				if(!hSim_){
					std::cerr << "Error retrieving sim histogram " << pixelhistoname << std::endl;
					std::cerr << "simFilePath = " << simFilePath << std::endl;
					return;
				}

				hSim_->SetDirectory(0);

				if((hData->GetNbinsX() != hSim_->GetNbinsX())){
					std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
					return;
				}

				dataIntegral = hData->Integral();
				simIntegral = hSim_->Integral();

				scaleFactor = 0.;

				if(simIntegral > 0){
					scaleFactor = dataIntegral/simIntegral;
				} else {
					std::cerr << "sim histogram has zero integral and thus cannot scale" << std::endl;
					return;
				}

				hSim_->Rebin(rebinfactor);
				hData->Rebin(rebinfactor);

				ndf = 0;
				chi2 = ComputeChi2BinByBin(hData, hSim_, scaleFactor, ndf);
				double pval = TMath::Prob(chi2,ndf);
				std::array<double,3> results = {chi2, ndf, pval};
				chi2_results.push_back(results);

				delete hSim_;
				hData->SetDirectory(0);
			}
			simfile->Close();

			double totalchi2 = 0.;
			double totalndf = 0.;
			for(const auto& res : chi2_results){
				totalchi2 += res[0];
				totalndf += res[1];
			}
			double totalreducedchi2 = totalchi2/totalndf;

			hBeamSpotSearch_4pixelspectra_chi2->Fill(xbin, ybin, totalchi2);
			hBeamSpotSearch_4pixelspectra_redchi2->Fill(xbin, ybin, totalreducedchi2);
		}
	}

	datafile->Close();
	//prep outfile:
	TString outfilename = "Lithium6_1plus_BeamSpotSearch.root";
	TFile *outfile = new TFile(outfilename,"RECREATE");

	//add labels for the hBeamSpotSearch_4pixel_redchi2 histograms
	int nx = hBeamSpotSearch_4pixel_redchi2->GetNbinsX();
	int ny = hBeamSpotSearch_4pixel_redchi2->GetNbinsY();

	std::vector<std::pair<double,int>> values;
	values.reserve(nx*ny);

	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){
			double val = hBeamSpotSearch_4pixel_redchi2->GetBinContent(ix,iy);
			int globalbin = hBeamSpotSearch_4pixel_redchi2->GetBin(ix,iy);
			if(val>0) values.push_back({val,globalbin});
		}
	}

	//sort by reduced chi2 value:
	std::sort(values.begin(), values.end(),
			  [](auto&a, auto&b){ return a.first < b.first; });

	std::map<int,int> rankmap;
	for(size_t i=0; i<values.size(); i++){
		rankmap[values[i].second] = i+1;
	}

	TLatex latex;
	latex.SetTextAlign(22);
	latex.SetTextSize(0.012);
	latex.SetTextColor(kBlack);

	TCanvas *c1 = new TCanvas("c1_4pix", "BeamspotSearch_4pix", 800, 800);
	c1->SetRightMargin(0.15);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);
	c1->SetTopMargin(0.08);

	hBeamSpotSearch_4pixel_redchi2->Draw("COLZ");
	hBeamSpotSearch_4pixel_redchi2->GetXaxis()->SetTitle("Beam #sigma_{x} (m)");
	hBeamSpotSearch_4pixel_redchi2->GetYaxis()->SetTitle("Beam #sigma_{y} (m)");
	hBeamSpotSearch_4pixel_redchi2->SetStats(0);

	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){

			double val = hBeamSpotSearch_4pixel_redchi2->GetBinContent(ix,iy);
			if(val <= 0) continue;

			int globalbin = hBeamSpotSearch_4pixel_redchi2->GetBin(ix,iy);
			int rank = rankmap[globalbin];

			double x = hBeamSpotSearch_4pixel_redchi2->GetXaxis()->GetBinCenter(ix);
			double y = hBeamSpotSearch_4pixel_redchi2->GetYaxis()->GetBinCenter(iy);

			double dy = 0.1*(hBeamSpotSearch_4pixel_redchi2->GetYaxis()->GetBinWidth(iy));

			 TString text = Form("#%d", rank);
			 latex.DrawLatex(x,y+dy,text);

			 text = Form("%.2f",val);
			 latex.DrawLatex(x,y-dy,text);

		}
	}




	gPad->SetFixedAspectRatio();
	outfile->cd();
	hBeamSpotSearch_4pixel_chi2->Write();
	hBeamSpotSearch_4pixel_redchi2->Write();
	c1->Write();
	c1->SaveAs("Lithium6_1plus_4pixel_redchi2.png");


	//add labels for the hBeamSpotSearch_16pixel_redchi2 histograms
	nx = hBeamSpotSearch_16pixel_redchi2->GetNbinsX();
	ny = hBeamSpotSearch_16pixel_redchi2->GetNbinsY();

	values.clear();
	values.reserve(nx*ny);

	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){
			double val = hBeamSpotSearch_16pixel_redchi2->GetBinContent(ix,iy);
			int globalbin = hBeamSpotSearch_16pixel_redchi2->GetBin(ix,iy);
			if(val>0) values.push_back({val,globalbin});
		}
	}

	//sort by reduced chi2 value:
	std::sort(values.begin(), values.end(),
			  [](auto&a, auto&b){ return a.first < b.first; });

	rankmap.clear();
	for(size_t i=0; i<values.size(); i++){
		rankmap[values[i].second] = i+1;
	}

	latex.SetTextAlign(22);
	latex.SetTextSize(0.012);
	latex.SetTextColor(kBlack);

	c1 = new TCanvas("c1_16pix", "BeamspotSearch_16pix", 800, 800);
	c1->SetRightMargin(0.15);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);
	c1->SetTopMargin(0.08);

	hBeamSpotSearch_16pixel_redchi2->Draw("COLZ");
	hBeamSpotSearch_16pixel_redchi2->GetXaxis()->SetTitle("Beam #sigma_{x} (m)");
	hBeamSpotSearch_16pixel_redchi2->GetYaxis()->SetTitle("Beam #sigma_{y} (m)");
	hBeamSpotSearch_16pixel_redchi2->SetStats(0);


	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){

			double val = hBeamSpotSearch_16pixel_redchi2->GetBinContent(ix,iy);
			if(val <= 0) continue;

			int globalbin = hBeamSpotSearch_16pixel_redchi2->GetBin(ix,iy);
			int rank = rankmap[globalbin];

			double x = hBeamSpotSearch_16pixel_redchi2->GetXaxis()->GetBinCenter(ix);
			double y = hBeamSpotSearch_16pixel_redchi2->GetYaxis()->GetBinCenter(iy);

			double dy = 0.1*(hBeamSpotSearch_16pixel_redchi2->GetYaxis()->GetBinWidth(iy));

			 TString text = Form("#%d", rank);
			 latex.DrawLatex(x,y+dy,text);

			 text = Form("%.2f",val);
			 latex.DrawLatex(x,y-dy,text);

		}
	}

	gPad->SetFixedAspectRatio();
	outfile->cd();
	hBeamSpotSearch_16pixel_chi2->Write();
	hBeamSpotSearch_16pixel_redchi2->Write();
	c1->Write();
	c1->SaveAs("Lithium6_1plus_16pixel_redchi2.png");




	//add labels for the hBeamSpotSearch_4pixelspectra_redchi2 histograms
	nx = hBeamSpotSearch_4pixelspectra_redchi2->GetNbinsX();
	ny = hBeamSpotSearch_4pixelspectra_redchi2->GetNbinsY();

	values.clear();
	values.reserve(nx*ny);

	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){
			double val = hBeamSpotSearch_4pixelspectra_redchi2->GetBinContent(ix,iy);
			int globalbin = hBeamSpotSearch_4pixelspectra_redchi2->GetBin(ix,iy);
			if(val>0) values.push_back({val,globalbin});
		}
	}

	//sort by reduced chi2 value:
	std::sort(values.begin(), values.end(),
			  [](auto&a, auto&b){ return a.first < b.first; });

	rankmap.clear();
	for(size_t i=0; i<values.size(); i++){
		rankmap[values[i].second] = i+1;
	}

	latex.SetTextAlign(22);
	latex.SetTextSize(0.012);
	latex.SetTextColor(kBlack);

	c1 = new TCanvas("c1_4pixenergy", "BeamspotSearch_4pixenergy", 800, 800);
	c1->SetRightMargin(0.15);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);
	c1->SetTopMargin(0.08);

	hBeamSpotSearch_4pixelspectra_redchi2->Draw("COLZ");
	hBeamSpotSearch_4pixelspectra_redchi2->GetXaxis()->SetTitle("Beam #sigma_{x} (m)");
	hBeamSpotSearch_4pixelspectra_redchi2->GetYaxis()->SetTitle("Beam #sigma_{y} (m)");
	hBeamSpotSearch_4pixelspectra_redchi2->SetStats(0);


	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){

			double val = hBeamSpotSearch_4pixelspectra_redchi2->GetBinContent(ix,iy);
			if(val <= 0) continue;

			int globalbin = hBeamSpotSearch_4pixelspectra_redchi2->GetBin(ix,iy);
			int rank = rankmap[globalbin];

			double x = hBeamSpotSearch_4pixelspectra_redchi2->GetXaxis()->GetBinCenter(ix);
			double y = hBeamSpotSearch_4pixelspectra_redchi2->GetYaxis()->GetBinCenter(iy);

			double dy = 0.1*(hBeamSpotSearch_4pixelspectra_redchi2->GetYaxis()->GetBinWidth(iy));

			 TString text = Form("#%d", rank);
			 latex.DrawLatex(x,y+dy,text);

			 text = Form("%.2f",val);
			 latex.DrawLatex(x,y-dy,text);

		}
	}

	gPad->SetFixedAspectRatio();
	outfile->cd();
	hBeamSpotSearch_4pixelspectra_chi2->Write();
	hBeamSpotSearch_4pixelspectra_redchi2->Write();
	c1->Write();
	c1->SaveAs("Lithium6_1plus_fourPixelEnergy_redchi2.png");

	

	outfile->Close();
	std::cout << "Outpout saved to " << outfilename << std::endl;

}

/*
---------------------------------------------------------------------------------------------------------
---------------								0 plus below 								  ---------------
---------------------------------------------------------------------------------------------------------
*/

TString Lithium6_0plus_sabrehits(TString beamstringx, TString beamstringy, TString anglestring, TString phistring, TString stragglestring){

	//TString anglestring = "16782278";

	//uncomment below for DESKTOP
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_0plus_output.root";
	// TString dataHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta%s_%sx_%sy_histos.root",anglestring.Data(),beamstringx.Data(),beamstringy.Data());
	// TString simHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	//uncomment below for LAPTOP
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_0plus_output.root";
	TString dataHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Li3563keV_7500keV_theta%s_phi_%s_%sx_%sy_%s_histos.root",anglestring.Data(), phistring.Data(), beamstringx.Data(),beamstringy.Data(),stragglestring.Data());
	TString simHistLocalPath = "SABRE/GEOM/hSABREARRAY_hitsMapLocal";

	TString pngoutname = Form("Lithium6_0plus_simVSdata_sabrehits_theta%s_phi_%s_%sx_%sy_%s.png",anglestring.Data(), phistring.Data(), beamstringx.Data(), beamstringy.Data(), stragglestring.Data());


	//prep output file name
	//TString outfile_root_name = Form("Lithium6_0plus_simVSdata_sabrehits_%s.root",beamstring.Data());


	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return "ERROR";
	}

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram" << std::endl;
		datafile->Close();
		return "ERROR";
	}
	hData->SetDirectory(0);
	datafile->Close();



	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return "ERROR";
	}

	TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
	if(!hSim){
		std::cerr << "Error retrieving sim histogram" << std::endl;
		simfile->Close();
		return "ERROR";
	}
	hSim->SetDirectory(0);
	simfile->Close();



	if((hData->GetNbinsX() != hSim->GetNbinsX()) || (hData->GetNbinsY() != hSim->GetNbinsY())){
		std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
		return "ERROR";
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
		return "ERROR";
	}

	//calculate chi2:
	double chi2 = hData->Chi2Test(hSim, "CHI2/NDF UW");


	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);



	hSim->Draw();

	TPaveText *pt = new TPaveText(0.65, 0.75, 0.88, 0.88, "NDC");
	pt->AddText(Form("%sx %sy (scale factor = %.3f)", beamstringx.Data(), beamstringy.Data(), scaleFactor));
	pt->AddText(Form("#chi^{2}/NDF = %.1f",chi2));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	TPaveText *ptangle = new TPaveText(0.12, 0.7, 0.27, 0.9, "NDC");
	TString temp;
	if(anglestring == "188208") temp = "19.8#circ #pm 1#circ";
	else if(anglestring == "178218") temp = "19.8#circ #pm 2#circ";
	else if(anglestring == "168228") temp = "19.8#circ #pm 3#circ";
	else if(anglestring == "158238") temp = "19.8#circ #pm 4#circ";
	else if(anglestring == "148248") temp = "19.8#circ #pm 5#circ";
	ptangle->AddText(Form("%s",temp.Data()));
	ptangle->AddText(Form("phi %s", phistring.Data()));
	ptangle->AddText(stragglestring);
	//ptangle->AddText(Form("#chi^{2}/NDF = %.1f",chi2));
	ptangle->SetFillColorAlpha(kWhite,0.5);
	ptangle->Draw();

	c1->Update();

	// TFile* outfile_root = new TFile(outfile_root_name,"RECREATE");
	// if(!outfile_root || outfile_root->IsZombie()){
	// 	std::cerr << "Error creating output file" << std::endl;
	// 	return;
	// }

	//hData->Write("hData_original");
	//hSim->Write("hSim_scaled");
	c1->SaveAs(pngoutname);
	c1->Clear();
	hData->Draw();
	pt = new TPaveText(0.65, 0.75, 0.88, 0.88, "NDC");
	pt->AddText("LiFha exp 1par ^{6}Li 1^{+} DATA");
	c1->Update();
	c1->SaveAs(Form("Lithium6_0plus_simVSdata_sabrehits_%s.png","1PAREXPDATA"));
	//outfile_root->Close();
	return pngoutname;
	//std::cout << "Finished! Output written to: " << outfile_root_name << std::endl;

}

void Lithium6_0plus_sabrehits_auto(){

	// std::vector<TString> beamstrings = {
	// 									"fixed",
	// 									"gaus001",
	// 									"gaus002",
	// 									"gaus003",
	// 									"gaus004",
	// 									"gaus005"
	// 									// "gaus006",
	// 									// "gaus007",
	// 									// "gaus008",
	// 									// "gaus009",
	// 									// "gaus010"
	// 								};

	std::vector<TString> beamstringsx = {
											"gaus001"
										};

	std::vector<TString> beamstringsy = {
											"gaus001"									
										};

	std::vector<TString> anglestrings = {
											//"193203",			// 19.8 +/- 0.5
											//"188208",			// 19.8 +/- 1.0
											"178218"			// 19.8 +/- 2.0
											//"168228",			// 19.8 +/- 3.0
											//"158238",			// 19.8 +/- 4.0
											//"148248"			// 19.8 +/- 5.0
										};

	TString straggle = "straggleOn";

	//std::vector<double> phis = {0.5, 1.0, 1.5, 2.0, 2.5};
	std::vector<double> phis = {1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0};
	std::vector<TString> phistrings;
	for(const auto& philower : phis){
		for(const auto& phiupper : phis){
			//phistrings.push_back(Form("%f_%f",-philower,phiupper));
			TString s = TString::Format("%s_%s", smartFormat(-philower).c_str(), smartFormat(phiupper).c_str());
			phistrings.push_back(s);
		}
	}



	for(const auto& bsx : beamstringsx){
		for(const auto& bsy : beamstringsy){
			std::vector<TString> pngs;
			for(const auto& angstr : anglestrings){
				for(const auto& phistr : phistrings){
					TString pngoutname = Lithium6_0plus_sabrehits(bsx, bsy, angstr, phistr, straggle);
					pngs.push_back(pngoutname);
				}
			}

			//for all anglestring sabre hit patterns for this given bsx, bsy combo let's make into a gif:
			TString gifname = Form("Lithium6_0plus_simVSdata_sabrehits_%sx_%sy.gif",bsx.Data(),bsy.Data());
			TString cmd = "convert -delay 100 -loop 0 ";
			for(const auto& fn : pngs){
				cmd += fn;
				cmd += " ";
			}
			cmd += gifname;

			std::cout << "Creating gif for bsx = " << bsx << " and bsy = " << bsy << std::endl;
			gSystem->Exec(cmd);
			std::cout << "Saved " << gifname << std::endl;

			//clean up pngs:
			cmd = "rm ";
			for(const auto& fn : pngs){
				cmd += fn;
				cmd += " ";
			}

			// gSystem->Exec(cmd);
		}
	}

}

TString Lithium6_0plus_pixelhistos(int ringChan, int wedgeChan, TString beamstringx, TString beamstringy, TString anglestring, TString phistring, TString stragglestring){

	int SABRE_ID = 3;

	TString pixelhistoname = Form("hSABRE%d_pixel_r%dw%d_ESummary",SABRE_ID,ringChan,wedgeChan);

	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_0plus_output.root";
	TString dataHistLocalPath = Form("SABRE/SABRE%d/Pixels/",SABRE_ID);
	dataHistLocalPath = dataHistLocalPath + pixelhistoname;

	TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Li3563keV_7500keV_theta%s_phi_%s_%sx_%sy_%s_histos.root",anglestring.Data(), phistring.Data(), beamstringx.Data(), beamstringy.Data(), stragglestring.Data());
	TString simHistLocalPath = Form("SABRE/SABRE%d/Pixels/",SABRE_ID);
	simHistLocalPath = simHistLocalPath + pixelhistoname;

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		std::cerr << "dataFilePath = " << dataFilePath << "\n\n";
		return "ERROR";
	}

	TH1 *hData = dynamic_cast<TH1*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram" << std::endl;
		std::cerr << "dataFilePath = " << dataFilePath << "\n";
		std::cerr << "dataHistLocalPath = " << dataHistLocalPath << "\n\n";
		datafile->Close();
		return "ERROR";
	}
	hData->SetDirectory(0);
	datafile->Close();

	TFile *simfile = new TFile(simFilePath,"READ");
	if(!simfile || simfile->IsZombie()){
		std::cerr << "Error opening sim file" << std::endl;
		return "ERROR";
	}

	TH1 *hSim = dynamic_cast<TH1*>(simfile->Get(simHistLocalPath));
	if(!hSim){
		std::cerr << "Error retrieving sim histogram" << std::endl;
		simfile->Close();
		return "ERROR";
	}
	hSim->SetDirectory(0);
	simfile->Close();

	if( hData->GetNbinsX() != hSim->GetNbinsX() ){
		std::cerr << "Histogram bin counts do not match on at least one axis!" << std::endl;
		return "ERROR";
	}

	double dataIntegral = hData->Integral();
	double simIntegral = hSim->Integral();

	double scaleFactor = 0.;

	if(simIntegral > 0){
		scaleFactor = dataIntegral/simIntegral;
		std::cout << "Global scale factor is " << scaleFactor << std::endl;
	} else {
		std::cerr << "sim histo has zero integral, cannot scale" << std::endl;
		return "ERROR";
	}

	const int rebinfactor = 4;

	double maxdata = hData->GetMaximum();
	double maxsim = hSim->GetMaximum();

	double ymax = std::max(maxdata,maxsim);

	double ndf;
	double chi2 = ComputeChi2BinByBin(hData, hSim, scaleFactor, ndf);

	TCanvas *c1 = new TCanvas("c1","Data vs Sim", 800, 600);
	gStyle->SetOptStat(0);

	hData->SetLineColor(kViolet);
	hData->SetLineWidth(4);
	hData->SetTitle(pixelhistoname);
	TH1D *hDataRebin = dynamic_cast<TH1D*>(hData->Rebin(rebinfactor));
	hDataRebin->GetXaxis()->SetRangeUser(0,3);
	hDataRebin->GetYaxis()->SetRangeUser(0,500);
	hDataRebin->Draw("HIST");

	hSim->SetLineColor(kOrange);
	hSim->SetLineWidth(2);
	TH1D* hSimRebin = dynamic_cast<TH1D*>(hSim->Rebin(rebinfactor));
	hSimRebin->GetXaxis()->SetRangeUser(0,3);
	hSimRebin->GetYaxis()->SetRangeUser(0,500);
	hSimRebin->Draw("HIST SAME");

	TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
	legend->AddEntry(hDataRebin,"Data","l");
	TString simlabel = Form("Sim (scaled x %.3f)",scaleFactor);
	legend->AddEntry(hSimRebin,simlabel,"l");
	legend->Draw();

	TPaveText *pt = new TPaveText(0.65, 0.63, 0.88, 0.74, "NDC");
	//pt->AddText(Form("%s",anglestring.Data()));
	pt->AddText(Form("%sx %sy", beamstringx.Data(), beamstringy.Data()));
	pt->AddText(Form("(scale factor = %.3f)", scaleFactor));
	pt->AddText(Form("#chi^{2}/NDF = %.1f",chi2/ndf));
	pt->SetFillColorAlpha(kWhite,0.5);
	pt->SetTextAlign(12);
	pt->Draw();

	TPaveText *ptangle = new TPaveText(0.12, 0.7, 0.27, 0.9, "NDC");
	TString temp;
	if(anglestring == "193203") temp = "19.8#circ #pm 0.5#circ";
	else if(anglestring == "188208") temp = "19.8#circ #pm 1#circ";
	else if(anglestring == "178218") temp = "19.8#circ #pm 2#circ";
	else if(anglestring == "168228") temp = "19.8#circ #pm 3#circ";
	else if(anglestring == "158238") temp = "19.8#circ #pm 4#circ";
	else if(anglestring == "148248") temp = "19.8#circ #pm 5#circ";
	ptangle->AddText(Form("%s",temp.Data()));
	//ptangle->AddText(Form("#chi^{2}/NDF = %.1f",chi2));
	ptangle->AddText(stragglestring);	
	ptangle->SetFillColorAlpha(kWhite,0.5);
	ptangle->Draw();

	c1->Update();

	TString pngoutname = Form("Lithium6_0plus_simVSdata_pixelhisto_r%dw%d_theta%s_%sx_%sy_%s.png",ringChan,wedgeChan,anglestring.Data(),beamstringx.Data(), beamstringy.Data(), stragglestring.Data());
	c1->SaveAs(pngoutname);

	delete c1;
	delete hData;
	//delete hDataRebin;
	delete hSim;
	//delete hSimRebin;

	return pngoutname;

}

void Lithium6_0plus_pixelhistos_auto(){

	TString beamx = "gaus0005";
	TString beamy = "gaus0005";

	TString anglestring = "178218";

	TString phistring = "-2.125_2.125";

	TString straggle = "straggleOff";

	std::vector<int> ringChans = {
									71,
									72,
									73	
								 };

	std::vector<int> wedgeChans = {
									29,
									30
								  };

	for(const auto& ring : ringChans){
		for(const auto& wedge : wedgeChans){

			Lithium6_0plus_pixelhistos(ring, wedge, beamx, beamy, anglestring, phistring, straggle);

		}
	}

	std::cout << "Done!" << std::endl;

}

void Lithium6_0plus_fourpixelchisquared(){

	TString anglestring = "178218";

	std::vector<TString> beamstrings = {
											"fixed"
											//"gaus001"
											// "gaus002",
											// "gaus003",
											// "gaus004",
											// "gaus005"
										};

	// TString bsx = "fixed";
	// TString bsy = "fixed";


	//std::vector<double> phis = {1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0};
	// std::vector<double> phis = {2.0, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3.0, 3.125, 3.25, 3.375, 3.5};
	// std::vector<TString> phistrings;
	// for(const auto& philower : phis){
	// 	for(const auto& phiupper : phis){
	// 		//phistrings.push_back(Form("%f_%f",-philower,phiupper));
	// 		TString s = TString::Format("%s_%s", smartFormat(-philower).c_str(), smartFormat(phiupper).c_str());
	// 		phistrings.push_back(s);
	// 	}
	// }

	std::vector<TString> phistrings = {"-3.5_3.0", "-2.125_2.125", "-3.0_2.75", "-3.5_2.625"};

	double xedges[14] = {-3.5625, -3.4375, -3.3125, -3.1875, -3.0625, -2.9375, -2.8125, -2.6875, -2.5625, -2.4375, -2.3125, -2.1875, -2.0625, -1.9375};
	double yedges[14] = {1.9375, 2.0625, 2.1875, 2.3125, 2.4375, 2.5625, 2.6875, 2.8125, 2.9375, 3.0625, 3.1875, 3.3125, 3.4375, 3.5625};

	TH2D *hGridSearchChi2 = new TH2D("hGridSearchChi2", "GridSearchChi2", 13, xedges, 13, yedges);
	TH2D *hGridSearchReducedChi2 = new TH2D("hGridSearchReducedChi2", "GridSearchReducedChi2", 13, xedges, 13, yedges);

	//---------------------------------------------------
	//				Establish data values
	//---------------------------------------------------

	//uncomment for DESKTOP:
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
	//uncomment for LAPTOP:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_0plus_output.root";
	TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";


	// TString path_pix_r71_w29 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r71_w29";
	// TString path_pix_r72_w29 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r72_w29";
	// TString path_pix_r72_w30 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r72_w30";
	// TString path_pix_r71_w30 = "SABRE/SABRE3/Pixels/hSABRE3_pixel_r71_w30";

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	// TH1 *hpix_r71_w29 = dynamic_cast<TH1*>(datafile->Get(path_pix_r71_w29));
	// TH1 *hpix_r72_w29 = dynamic_cast<TH1*>(datafile->Get(path_pix_r72_w29));
	// TH1 *hpix_r72_w30 = dynamic_cast<TH1*>(datafile->Get(path_pix_r72_w30));
	// TH1 *hpix_r71_w30 = dynamic_cast<TH1*>(datafile->Get(path_pix_r71_w30));
	// if(!hpix_r71_w29 || !hpix_r72_w29 || !hpix_r72_w30 || !hpix_r71_w30){
	// 	std::cerr << "Error retrieving at least one data histogram!" << std::endl;
	// 	return;
	// }

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram!" << std::endl;
		return;
	}

	// hpix_r71_w29->SetDirectory(0);
	// hpix_r72_w29->SetDirectory(0);
	// hpix_r72_w30->SetDirectory(0);
	// hpix_r71_w30->SetDirectory(0);
	hData->SetDirectory(0);
	datafile->Close();

	// double counts_pix_r71_w29 = hpix_r71_w29->GetEntries();
	// double counts_pix_r72_w29 = hpix_r72_w29->GetEntries();
	// double counts_pix_r72_w30 = hpix_r72_w30->GetEntries();
	// double counts_pix_r71_w30 = hpix_r71_w30->GetEntries();

	double dataIntegral = (hData->GetBinContent(hData->GetBin(30,8))) + (hData->GetBinContent(hData->GetBin(30,9))) + (hData->GetBinContent(hData->GetBin(31,9))) + (hData->GetBinContent(hData->GetBin(31,8)));

	//std::vector<double> fourpix_relcounts = {counts_pix_r71_w29/fourpixsum, counts_pix_r72_w29/fourpixsum, counts_pix_r72_w30/fourpixsum, counts_pix_r71_w30/fourpixsum};

	//for(const auto& fn : filenames){

	for(const auto& ps : phistrings){
		for(const auto& bsx : beamstrings){
			for(const auto& bsy : beamstrings){

				TString fn = Form("kin2mc_7Li3He4He6Li3563keV_7500keV_theta%s_phi_%s_%sx_%sy_histos.root", anglestring.Data(), ps.Data(), bsx.Data(), bsy.Data());

				//establish sim values:
				//uncomment for DESKTOP
				// TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/%s",fn.Data());
				// TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
				//uncomment for LAPTOP:
				TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/%s",fn.Data());
				TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";


				TFile *simfile = new TFile(simFilePath,"READ");
				if(!simfile || simfile->IsZombie()){
					std::cerr << "Error opening sim file " << simFilePath << "\n\n";
					continue;
				}


				TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
				if(!hSim){
					std::cerr << "Error retrieving sim histo from " << simFilePath << "\n\n";
					continue;
				}

				double simIntegral = (hSim->GetBinContent(hSim->GetBin(30,8))) + (hSim->GetBinContent(hSim->GetBin(30,9))) + (hSim->GetBinContent(hSim->GetBin(31,9))) + (hSim->GetBinContent(hSim->GetBin(31,8)));


				hSim->SetDirectory(0);
				simfile->Close();

				if( hSim->GetNbinsX() != hData->GetNbinsX() ){
					std::cerr << "Histogram binning does not match for data/sim histogram pair\n";
					continue;
				}




				double scaleFactor;
				if(simIntegral > 0){
					scaleFactor = dataIntegral/simIntegral;
					hSim->Scale(scaleFactor);
				} else {
					std::cerr << "sim histogram from " << fn.Data() << " has zero integral and thus cannot be scaled! Continuing..." << std::endl;
					continue;
				}

				//calculate the chi squared for this sim file:

				double r73w29_data = (hData->GetBinContent(hData->GetBin(30,10)));//binx = 30 is wedge 29, biny = 10 is local ring 9 which is ring 73 for SABRE3
				double r72w29_data = (hData->GetBinContent(hData->GetBin(30,9)));//binx = 30 is wedge 29, biny = 9 is local ring 8 which is ring 72 for SABRE3
				double r73w30_data = (hData->GetBinContent(hData->GetBin(31,10)));//binx = 31 is wedge 30, biny = 10 is local ring 9 which is ring 73 for SABRE3
				double r72w30_data = (hData->GetBinContent(hData->GetBin(31,9)));//binx = 31 is wedge 30, biny = 9 is local ring 8 for which is ring 72 for SABRE3

				double r73w29_data_error = (hData->GetBinError(hData->GetBin(30,10)));
				double r72w29_data_error = (hData->GetBinError(hData->GetBin(30,9)));
				double r73w30_data_error = (hData->GetBinError(hData->GetBin(31,10)));
				double r72w30_data_error = (hData->GetBinError(hData->GetBin(31,9)));

				double r73w29_sim = (hSim->GetBinContent(hSim->GetBin(30,10)));
				double r72w29_sim = (hSim->GetBinContent(hSim->GetBin(30,9)));
				double r73w30_sim = (hSim->GetBinContent(hSim->GetBin(31,10)));
				double r72w30_sim = (hSim->GetBinContent(hSim->GetBin(31,9)));

				double r73w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,10)));
				double r72w29_sim_error = (hSim->GetBinError(hSim->GetBin(30,9)));
				double r73w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,10)));
				double r72w30_sim_error = (hSim->GetBinError(hSim->GetBin(31,9)));

				double chi2 = ( std::pow((r73w29_data - r73w29_sim),2) / ( std::pow(r73w29_data_error,2) + std::pow(r73w29_sim_error,2)) ) + 
							  ( std::pow((r72w29_data - r72w29_sim),2) / ( std::pow(r72w29_data_error,2) + std::pow(r72w29_sim_error,2)) ) +
							  ( std::pow((r73w30_data - r73w30_sim),2) / ( std::pow(r73w30_data_error,2) + std::pow(r73w30_sim_error,2)) ) +
							  ( std::pow((r72w30_data - r72w30_sim),2) / ( std::pow(r72w30_data_error,2) + std::pow(r72w30_sim_error,2)) );

				Int_t ndf = 4 - 1;//4 pixels minus the 1 DOF (the scale factor is only fitted parameter here)

				double reducedchi2 = chi2/ndf;


				//fill chi2 histogram from grid search:
				int xbin = 0;
				int ybin = 0;

				// if(bsx == "fixed") xbin = 0;
				// else if(bsx == "gaus001") xbin = 1;
				// else if(bsx == "gaus002") xbin = 2;
				// else if(bsx == "gaus003") xbin = 3;
				// else if(bsx == "gaus004") xbin = 4;
				// else if(bsx == "gaus005") xbin = 5;

				// if(bsy == "fixed") ybin = 0;
				// else if(bsy == "gaus001") ybin = 1;
				// else if(bsy == "gaus002") ybin = 2;
				// else if(bsy == "gaus003") ybin = 3;
				// else if(bsy == "gaus004") ybin = 4;
				// else if(bsy == "gaus005") ybin = 5;

				std::pair<TString,TString> philowhigh = parsephistring(ps);
				TString philow = philowhigh.first;
				TString phihigh = philowhigh.second;

				// if(philow == "-2.5") xbin = 5;
				// else if(philow == "-2.0") xbin = 4;
				// else if(philow == "-1.5") xbin = 3;
				// else if(philow == "-1.0") xbin = 2;
				// else if(philow == "-0.5") xbin = 1;

				// if(phihigh == "2.5") ybin = 5;
				// else if(phihigh == "2.0") ybin = 4;
				// else if(phihigh == "1.5") ybin = 3;
				// else if(phihigh == "1.0") ybin = 2;
				// else if(phihigh == "0.5") ybin = 1;

				hGridSearchChi2->Fill(philow.Atof(),phihigh.Atof(),chi2);
				hGridSearchReducedChi2->Fill(philow.Atof(),phihigh.Atof(),reducedchi2);
			}
		}
	}

	//add labels showing value and relative rank on each bin:
	int nx = hGridSearchReducedChi2->GetNbinsX();
	int ny = hGridSearchReducedChi2->GetNbinsY();

	//store (value, global_bin) for ranking
	std::vector<std::pair<double,int>> values;
	values.reserve(nx*ny);

	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){
			double val = hGridSearchReducedChi2->GetBinContent(ix,iy);
			int globalbin = hGridSearchReducedChi2->GetBin(ix,iy);

			if(val > 0) values.push_back({val,globalbin});
		}
	}

	//sort by reduced chi2 value:
	std::sort(values.begin(), values.end(),
		      [](auto&a, auto&b){ return a.first < b.first; });

	std::map<int,int> rankmap;
	for(size_t i=0; i<values.size(); i++){
		rankmap[values[i].second] = i+1;
	}

	TLatex latex;
	latex.SetTextAlign(22);
	latex.SetTextSize(0.012);
	latex.SetTextColor(kBlack);

	TFile *outfile = new TFile("GridSearch.root","RECREATE");

	outfile->cd();
	//hData->Write();
	//hSim->Write();
	hGridSearchChi2->Write();
	hGridSearchReducedChi2->Write();

	//set up square TCanvas here:
	TCanvas *c1 = new TCanvas("c1","Grid Search Reduced Chi2",800,800);
	c1->SetRightMargin(0.15);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);
	c1->SetTopMargin(0.08);

	hGridSearchReducedChi2->Draw("COLZ");
	hGridSearchReducedChi2->GetXaxis()->SetTitle("#phi_{low} (#circ)");
	hGridSearchReducedChi2->GetYaxis()->SetTitle("#phi_{high} (#circ)");
	hGridSearchReducedChi2->SetStats(0);

	//draw labels on each bin
	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){

			double val = hGridSearchReducedChi2->GetBinContent(ix,iy);
			if(val <= 0) continue;

			int globalbin = hGridSearchReducedChi2->GetBin(ix,iy);
			int rank = rankmap[globalbin];

			double x = hGridSearchReducedChi2->GetXaxis()->GetBinCenter(ix);
			double y = hGridSearchReducedChi2->GetYaxis()->GetBinCenter(iy);

			double dy = 0.1*(hGridSearchReducedChi2->GetYaxis()->GetBinWidth(iy));

			TString text = Form("#%d",rank);
			latex.DrawLatex(x,y+dy,text);

			text = Form("%.2f",val);
			latex.DrawLatex(x,y-dy,text);

		}
	}

	//force equal scaling on both x and y axis w/ gpad:
	gPad->SetFixedAspectRatio();
	outfile->cd();
	c1->Write("cGridSearchReducedChi2");
	c1->SaveAs("Lithium6_0plus_fourpixelchisquared.png");

	outfile->Close();

}

void Lithium6_0plus_sixteenpixelchisquared(){

	TString anglestring = "178218";

	std::vector<TString> beamstrings = {
											//"fixed"
											"gaus001"
											// "gaus002",
											// "gaus003",
											// "gaus004",
											// "gaus005"
										};

	//std::vector<double> phis = {1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0};
	std::vector<double> phis = {2.0, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3.0, 3.125, 3.25, 3.375, 3.5};								
	std::vector<TString> phistrings;
	for(const auto& philower : phis){
		for(const auto& phiupper : phis){
			//phistrings.push_back(Form("%f_%f",-philower,phiupper));
			TString s = TString::Format("%s_%s", smartFormat(-philower).c_str(), smartFormat(phiupper).c_str());
			phistrings.push_back(s);
		}
	}

	// double xedges[8] = {-3.125, -2.875, -2.625, -2.375, -2.125, -1.875, -1.625, -1.375};
	// double yedges[8] = {1.375, 1.625, 1.875, 2.125, 2.375, 2.625, 2.875, 3.125};

	double xedges[14] = {-3.5625, -3.4375, -3.3125, -3.1875, -3.0625, -2.9375, -2.8125, -2.6875, -2.5625, -2.4375, -2.3125, -2.1875, -2.0625, -1.9375};
	double yedges[14] = {1.9375, 2.0625, 2.1875, 2.3125, 2.4375, 2.5625, 2.6875, 2.8125, 2.9375, 3.0625, 3.1875, 3.3125, 3.4375, 3.5625};

	TH2D *hGridSearchChi2 = new TH2D("hGridSearchChi2", "GridSearchChi2", 13, xedges, 13, yedges);
	TH2D *hGridSearchReducedChi2 = new TH2D("hGridSearchReducedChi2", "GridSearchReducedChi2", 13, xedges, 13, yedges);


	//---------------------------------------------------
	//				Establish data values
	//---------------------------------------------------

	//uncomment for DESKTOP:
	// TString dataFilePath = "/home/zmpur/SABREsim/det/ROOT/LiFha_1par_exp_1plus_output.root";
	// TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";
	//uncomment for LAPTOP:
	TString dataFilePath = "/mnt/e/RMSRecon/etc/zmpROOT/LiFha_1par_exp_0plus_output.root";
	TString dataHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

	TFile *datafile = new TFile(dataFilePath,"READ");
	if(!datafile || datafile->IsZombie()){
		std::cerr << "Error opening data file" << std::endl;
		return;
	}

	TH2 *hData = dynamic_cast<TH2*>(datafile->Get(dataHistLocalPath));
	if(!hData){
		std::cerr << "Error retrieving data histogram!" << std::endl;
		return;
	}

	hData->SetDirectory(0);
	datafile->Close();

	//ring (y) range for 3.563 MeV 0+ state from data
	//same with wedge/x range but this does not change between 1+ GS and 0+ 3.563 MeV state
	int x_minbin = 28+1;//add 1 because root bin 1 is corresponding to our ring/wedge 0, there fore our ring/wedge N is bin N+1
	int x_maxbin = 31+1;//add 1 because root bin 1 is corresponding to our ring/wedge 0, there fore our ring/wedge N is bin N+1
	int y_minbin = 7+1;//add 1 because root bin 1 is corresponding to our ring/wedge 0, there fore our ring/wedge N is bin N+1
	int y_maxbin = 10+1;//add 1 because root bin 1 is corresponding to our ring/wedge 0, there fore our ring/wedge N is bin N+1
	double dataIntegral = hData->Integral(x_minbin, x_maxbin, y_minbin, y_maxbin);

	for(const auto& ps : phistrings){
		for(const auto& bsx : beamstrings){
			for(const auto& bsy : beamstrings){

				TString fn = Form("kin2mc_7Li3He4He6Li3563keV_7500keV_theta%s_phi_%s_%sx_%sy_histos.root", anglestring.Data(), ps.Data(), bsx.Data(), bsy.Data());

				//uncomment for DESKTOP:
				//TString simFilePath = Form("/home/zmpur/SABREsim/det/kin2mc/%s",fn.Data());
				//TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

				//uncomment for LAPTOP:
				TString simFilePath = Form("/mnt/e/SABREsim/det/kin2mc/%s",fn.Data());
				TString simHistLocalPath = "SABRE/hSABRE_RingsVSWedges";

				TFile *simfile = new TFile(simFilePath,"READ");
				if(!simfile || simfile->IsZombie()){
					std::cerr << "Error opening sim file " << simFilePath << "\n\n";
					continue;
				}

				TH2 *hSim = dynamic_cast<TH2*>(simfile->Get(simHistLocalPath));
				if(!hSim){
					std::cerr << "Error retrieving sim histo from " << simFilePath << "\n\n";
					continue;
				}

				double simIntegral = hSim->Integral(x_minbin, x_maxbin, y_minbin, y_maxbin);

				hSim->SetDirectory(0);
				simfile->Close();

				if(hSim->GetNbinsX() != hData->GetNbinsX()){
					std::cerr << "Histogram binning does not match for data/sim histogram pair\n\n";
					continue;
				}

				double scaleFactor;
				if(simIntegral > 0){
					scaleFactor = dataIntegral/simIntegral;
					hSim->Scale(scaleFactor);
				} else {
					std::cerr << "sim histo from " << fn.Data() << " has zero integral and thus cannot be scaled!\n\n";
					continue;
				}

				double chi2 = 0.;
				double reducedchi2 = 0.;
				int ndf = 16 - 1;//16 pixels minus 1 DOF (scale factor only DOF on a per-run basis)

				for(int x=x_minbin; x<=x_maxbin; x++){
					
					for(int y=y_minbin; y<=y_maxbin; y++){

						double datavalue = hData->GetBinContent(hData->GetBin(x,y));
						double simvalue = hSim->GetBinContent(hSim->GetBin(x,y));

						if(datavalue == 0 && simvalue == 0){
							ndf -= 1;
							std::cout << "passing for x = " << x << " and y = " << y << "\n";
							continue;
						}
						
						double dataerror = hData->GetBinError(hData->GetBin(x,y));
						double simerror = hSim->GetBinError(hSim->GetBin(x,y));

						double numerator = std::pow((datavalue - simvalue) , 2);
						double denominator = std::pow(dataerror, 2) + std::pow(simerror, 2);
						chi2 += (numerator/denominator);

					}
				}

				reducedchi2 = chi2/ndf;

				int xbin = 0;
				int ybin = 0;

				// if(bsx == "fixed") xbin = 0;
				// else if(bsx == "gaus001") xbin = 1;
				// else if(bsx == "gaus002") xbin = 2;
				// else if(bsx == "gaus003") xbin = 3;
				// else if(bsx == "gaus004") xbin = 4;
				// else if(bsx == "gaus005") xbin = 5;

				// if(bsy == "fixed") ybin = 0;
				// else if(bsy == "gaus001") ybin = 1;
				// else if(bsy == "gaus002") ybin = 2;
				// else if(bsy == "gaus003") ybin = 3;
				// else if(bsy == "gaus004") ybin = 4;
				// else if(bsy == "gaus005") ybin = 5;

				std::pair<TString,TString> philowhigh = parsephistring(ps);
				TString philow = philowhigh.first;
				TString phihigh = philowhigh.second;

				hGridSearchChi2->Fill(philow.Atof(),phihigh.Atof(),chi2);
				hGridSearchReducedChi2->Fill(philow.Atof(),phihigh.Atof(),reducedchi2);
			}
		}
	}

	//add labels showing values and relative ranks on each bin:
	int nx = hGridSearchReducedChi2->GetNbinsX();
	int ny = hGridSearchReducedChi2->GetNbinsY();

	//store (value, global_bin) for ranking
	std::vector<std::pair<double,int>> values;
	values.reserve(nx*ny);

	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){
			double val = hGridSearchReducedChi2->GetBinContent(ix,iy);
			int globalbin = hGridSearchReducedChi2->GetBin(ix,iy);

			if(val > 0) values.push_back({val,globalbin});
		}
	}

	//sort by reduced chi2 value:
	std::sort(values.begin(), values.end(),
		      [](auto&a, auto&b){ return a.first < b.first; });

	std::map<int,int> rankmap;
	for(size_t i=0; i<values.size(); i++){
		rankmap[values[i].second] = i+1;
	}

	TLatex latex;
	latex.SetTextAlign(22);
	latex.SetTextSize(0.012);
	latex.SetTextColor(kBlack);

	TFile *outfile = new TFile("GridSearch_16pix.root","RECREATE");
	outfile->cd();
	hGridSearchChi2->Write();
	hGridSearchReducedChi2->Write();

	TCanvas *c1 = new TCanvas("c1","Grid Seach Reduced Chi2",800,800);
	c1->SetRightMargin(0.15);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);
	c1->SetTopMargin(0.08);

	hGridSearchReducedChi2->Draw("COLZ");
	hGridSearchReducedChi2->GetXaxis()->SetTitle("#phi_{low} (#circ)");
	hGridSearchReducedChi2->GetYaxis()->SetTitle("#phi_{high} (#circ)");
	hGridSearchReducedChi2->SetStats(0);

	//draw labels on each bin:
	for(int ix=1; ix<=nx; ix++){
		for(int iy=1; iy<=ny; iy++){

			double val = hGridSearchReducedChi2->GetBinContent(ix,iy);
			if(val <= 0) continue;

			int globalbin = hGridSearchReducedChi2->GetBin(ix,iy);
			int rank = rankmap[globalbin];

			double x = hGridSearchReducedChi2->GetXaxis()->GetBinCenter(ix);
			double y = hGridSearchReducedChi2->GetYaxis()->GetBinCenter(iy);

			double dy = 0.1*(hGridSearchReducedChi2->GetYaxis()->GetBinWidth(iy));

			TString text = Form("#%d",rank);
			latex.DrawLatex(x,y+dy,text);

			text = Form("%.2f",val);
			latex.DrawLatex(x,y-dy,text);

		}
	}

	gPad->SetFixedAspectRatio();
	outfile->cd();
	c1->Write("cGridSearchReducedChi2");
	c1->SaveAs("Lithium6_0plus_sixteenpixelchisquared.png");

	outfile->Close();
}

/*
---------------------------------------------------------------------------------------------------------
---------------								3 plus below 								  ---------------
---------------------------------------------------------------------------------------------------------
*/

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


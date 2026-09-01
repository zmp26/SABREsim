#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TH1D.h>
#include <TDirectory.h>

//identical to EffCor_9B_p8Be except for effFile default changed (and exact mapping for effIndex slightly different, since different decay topology)
void EffCor_9B_a5Li(const char* histName = "h", int effIndex = 9, const char* effFile = "/home/zachpurcell/efficiencies/9B_4He_4He_eff_1690_6000_SABREres50keV_thresh100keV.txt"){

	//histName -> name in ROOT session of histogram to be efficiency-corrected (default is "h", but can be whatever name you choose)
	//effFile  -> absolute path to txt file containing efficiency parameters for 9B->a1+5Li->p+a2
	//effIndex -> following are the choices for the particular efficiency (column) in the effFile:
	//		 EffIndex | Particles Detected						  |		Meaning of Efficiency
	//			1 	  =		1 										=====>	efficiency for detecting first decay only (a1 for this case)
	//			2 	  =		2 										=====>	efficiency for detecting second decay only (p for this case)
	//			3 	  =		3 										=====>	efficiency for detecting third decay only (a2 for this case)
	//			4 	  = 	any 1 particle 							=====>	efficiency for detecting any 1 particle
	//			5 	  =		12 || 21 								=====>	efficiency for detecting first decay and second (a1, p for this case)
	//			6 	  =		23 || 32 								=====>	efficiency for detecting second decay and last (p, a2 for this case)
	//			7 	  =		13 || 31 								=====>	efficiency for detecting first decay and last (a1, a2 for this case)
	//			8	  =		any 2 particles							=====>	efficiency for detecting any two particles
	//			9	  =		123 || 132 || 213 || 231 || 312 || 321	=====>	efficiency for detecting all three resonance decay particles
	//			
	//		default for effIndex is 9 for no particular reason, so make sure to change as needed during call

	if(effIndex < 1 || effIndex > 9){
		std::cout << "Expect effIndex integer in range 1-9" << std::endl;
		return;
	}


	TH1D *h_uncorr = dynamic_cast<TH1D*>(gDirectory->Get(histName));
	if(!h_uncorr){
		std::cerr << "Error: histogram with name " << histName << " not found in current session!" << std::endl;
		return;
	}

	int nbins = h_uncorr->GetNbinsX();
	double xmin = h_uncorr->GetXaxis()->GetXmin();
	double xmax = h_uncorr->GetXaxis()->GetXmax();
	double binWidth_keV = h_uncorr->GetBinWidth(1)*1000.;

	if (std::abs(binWidth_keV - 5.0) > 0.01) {
		std::cout << "Warning: Target histogram bin width is " << binWidth_keV << " keV/bin (efficiency file assumes 5.0 keV/bin)." << std::endl;
	}


	TH1D *h_eff = new TH1D(Form("%s_eff_temp", histName), "Efficiency Histo", nbins, xmin, xmax);

	std::ifstream infile(effFile);
	if(!infile.is_open()){
		std::cerr << "Error: Could not open efficiency file " << effFile << "!" << std::endl;
		delete h_eff;
		return;
	}

	std::string line;
	std::getline(infile, line);//skips header line

	double energy_keV;
	while(std::getline(infile, line)){
		if(line.empty()) continue;
		std::stringstream ss(line);
		if(!(ss >> energy_keV)) continue;

		double bu1=0., bu2=0., bu3=0., bu1bu2=0., bu2bu3=0., bu1bu3=0., bu1bu2bu3=0.;
		if(ss >> bu1 >> bu2 >> bu3 >> bu1bu2 >> bu2bu3 >> bu1bu3 >> bu1bu2bu3){
			double eff_percent = 0.;

			switch(effIndex){
				case 1: eff_percent = bu1; break;
				case 2: eff_percent = bu2; break;
				case 3: eff_percent = bu3; break;
				case 4: eff_percent = bu1 + bu2 + bu3; break;
				case 5: eff_percent = bu1bu2; break;
				case 6: eff_percent = bu2bu3; break;
				case 7: eff_percent = bu1bu3; break;
				case 8: eff_percent = bu1bu2 + bu2bu3 + bu1bu3; break;
				case 9: eff_percent = bu1bu2bu3; break;
			}

			double energy_MeV = energy_keV / 1000.;
			double eff_fraction = eff_percent/100.;

			int bin = h_eff->FindBin(energy_MeV + 0.5*h_eff->GetBinWidth(1));

			if(bin >= 1 && bin <= nbins){
				h_eff->SetBinContent(bin, eff_fraction);
				h_eff->SetBinError(bin, 0.);
			}
		}
	}

	infile.close();

	TString corrName = Form("%s_corr", histName);
	TH1D *h_corr = dynamic_cast<TH1D*>(h_uncorr->Clone(corrName));
	h_corr->SetTitle(Form("%s (Efficiency Corrected)", h_uncorr->GetTitle()));
	h_corr->Sumw2();

	//prevent divide by 0 errors just in case
	int zeroEffBins = 0;
	for(int i=1; i<=nbins; i++){
		double uncorr_val = h_uncorr->GetBinContent(i);
		double uncorr_err = h_uncorr->GetBinError(i);
		double eff_val = h_eff->GetBinContent(i);

		if(eff_val > 0.){
			h_corr->SetBinContent(i, uncorr_val/eff_val);
			h_corr->SetBinError(i, uncorr_err/eff_val);
		} else {
			h_corr->SetBinContent(i, 0.);
			h_corr->SetBinError(i, 0.);
			if(uncorr_val > 0.) zeroEffBins += 1;
		}
	}

	delete h_eff;

	if(zeroEffBins > 0){
		std::cout << "Warning: " << zeroEffBins << " bin(s) had counts with efficiency of 0%. These bins were set to 0." << std::endl;
	}

	h_corr->SetLineColor(kRed);
	h_corr->Draw("HIST SAME");

	std::cout << "Successfully created efficiency-corrected histogram: '" << corrName << "'" << std::endl;
}

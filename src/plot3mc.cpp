#include "plot3mc.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Vec3.h"

plot3mc::plot3mc(const std::string& outname){

	rand_gen = new TRandom3(0);

	sabre_thetaphimap = readAngleMaps();

	histofile = new TFile(outname.c_str(), "RECREATE");
	histoman = new HistoManager(histofile);
	histoman->loadHistoConfig("config/_3body.HMConfig");

	masstable = new MassTable;
	masstable->Init("config/masstable.dat");

}

plot3mc::~plot3mc(){
	delete histofile;
	delete histoman;
	delete masstable;
	delete rand_gen;
}

void plot3mc::FillKinematicsHistos(PHYSDATA& pd1, PHYSDATA& pd2, PHYSDATA& pd3, PHYSDATA& pd4){
	//basic kinematics histograms:
	histoman->getHisto1D("hELab1")->Fill(pd1.e);
	histoman->getHisto1D("hELab2")->Fill(pd2.e);
	histoman->getHisto1D("hELab3")->Fill(pd3.e);
	histoman->getHisto1D("hELab4")->Fill(pd4.e);
	histoman->getHisto1D("hThetaLab1")->Fill(pd1.theta);
	histoman->getHisto1D("hThetaLab2")->Fill(pd2.theta);
	histoman->getHisto1D("hThetaLab3")->Fill(pd3.theta);
	histoman->getHisto1D("hThetaLab4")->Fill(pd4.theta);
	histoman->getHisto1D("hPhiLab1")->Fill(pd1.phi);
	histoman->getHisto1D("hPhiLab2")->Fill(pd2.phi);
	histoman->getHisto1D("hPhiLab3")->Fill(pd3.phi);
	histoman->getHisto1D("hPhiLab4")->Fill(pd4.phi);
	histoman->getHisto2D("hELabThetaLab_1")->Fill(pd1.theta,pd1.e);
	histoman->getHisto2D("hELabThetaLab_2")->Fill(pd2.theta,pd2.e);
	histoman->getHisto2D("hELabThetaLab_3")->Fill(pd3.theta,pd3.e);
	histoman->getHisto2D("hELabThetaLab_4")->Fill(pd4.theta,pd4.e);
	histoman->getHisto2D("hELabPhiLab_1")->Fill(pd1.phi,pd1.e);
	histoman->getHisto2D("hELabPhiLab_2")->Fill(pd2.phi,pd2.e);
	histoman->getHisto2D("hELabPhiLab_3")->Fill(pd3.phi,pd3.e);
	histoman->getHisto2D("hELabPhiLab_4")->Fill(pd4.phi,pd4.e);
	histoman->getHisto2D("hThetaLabPhiLab_1")->Fill(pd1.theta,pd1.phi);
	histoman->getHisto2D("hThetaLabPhiLab_2")->Fill(pd2.theta,pd2.phi);
	histoman->getHisto2D("hThetaLabPhiLab_3")->Fill(pd3.theta,pd3.phi);
	histoman->getHisto2D("hThetaLabPhiLab_4")->Fill(pd4.theta,pd4.phi);

	TVector3 lab1(sin(DEGRAD*pd3.theta)*cos(DEGRAD*pd3.phi), sin(DEGRAD*pd3.theta)*sin(DEGRAD*pd3.phi), cos(DEGRAD*pd3.theta)), lab2(sin(DEGRAD*pd4.theta)*cos(DEGRAD*pd4.phi), sin(DEGRAD*pd4.theta)*sin(DEGRAD*pd4.phi), cos(DEGRAD*pd4.theta));
	histoman->getHisto1D("hBreakupLabAngle")->Fill(lab1.Angle(lab2)/DEGRAD);
}

void plot3mc::FillBeamSpotHisto(Vec3& reactionOrigin){
	histoman->getHisto2D("hBeamSpot")->Fill(reactionOrigin.GetX(),reactionOrigin.GetY());
}

void plot3mc::FillSABREHistos(SABREDATA& sd1, PHYSDATA& pd1){

	int ringoffset = offsets[sd1.detectorIndex].first;
	int wedgeoffset = offsets[sd1.detectorIndex].second;
	int globalring = sd1.ring + ringoffset;
	int globalwedge = sd1.wedge + wedgeoffset;

	TString pixelhistoname = Form("hSABRE%d_pixel_r%dw%d_ESummary", sd1.detectorIndex, globalring, globalwedge);
	histoman->getHisto1D(pixelhistoname)->Fill(sd1.ringEnergy);

	//rings vs wedges histogram:
	int zeroToFortyWedge = sd1.detectorIndex*numwedges + sd1.wedge;//0-7 for SABRE0, 8-15 for SABRE1, etc.
	histoman->getHisto2D("hSABRE_RingsVSWedges")->Fill(zeroToFortyWedge,sd1.ring);

	//SABRE vs kin2mc
	histoman->getHisto2D("h_SABRETheta_vs_kinTheta")->Fill(pd1.theta,sd1.theta);
	histoman->getHisto2D("h_SABREPhi_vs_kinPhi")->Fill(pd1.phi,sd1.phi);
	histoman->getHisto2D("h_SABRERingE_vs_kinE")->Fill(pd1.e,sd1.ringEnergy);
	histoman->getHisto2D("h_SABREWedgeE_vs_kinE")->Fill(pd1.e,sd1.wedgeEnergy);

	histoman->getHisto2D("hSABREARRAY_hitsMapLocal")->Fill(-sd1.localx, -sd1.localy);

	//SABRE histos
	TString hname = Form("hSABRE%d_RingHit",sd1.detectorIndex);
	histoman->getHisto1D(hname)->Fill(sd1.ring);


	hname = Form("hSABRE%d_WedgeHit",sd1.detectorIndex);
	histoman->getHisto1D(hname)->Fill(sd1.wedge);


	hname = Form("hSABRE%d_ESummaryWedges",sd1.detectorIndex);
	histoman->getHisto2D(hname)->Fill(sd1.wedge, sd1.wedgeEnergy);


	hname = Form("hSABRE%d_ESummaryRings",sd1.detectorIndex);
	histoman->getHisto2D(hname)->Fill(sd1.ring, sd1.ringEnergy);


	hname = Form("hSABRE%d_hitsMapLocal",sd1.detectorIndex);
	histoman->getHisto2D(hname)->Fill(-sd1.localx, -sd1.localy);


	hname = Form("hSABRE_AngleHitsMap");
	histoman->getHisto2D(hname)->Fill(sd1.theta, sd1.phi);


	hname = Form("hSABREARRAY_hitsMapLocal");
	histoman->getHisto2D(hname)->Fill(-sd1.localx, -sd1.localy);


	hname = Form("hSABRE%d_ringEVSkinE",sd1.detectorIndex);
	histoman->getHisto2D(hname)->Fill(pd1.e, sd1.ringEnergy);


	hname = Form("hSABRE%d_wedgeEVSkinE",sd1.detectorIndex);
	histoman->getHisto2D(hname)->Fill(pd1.e, sd1.wedgeEnergy);


	hname = Form("hSABRE%d_ringThetaVSkinTheta",sd1.detectorIndex);
	histoman->getHisto2D(hname)->Fill(pd1.theta, sd1.theta);


	hname = Form("hSABRE%d_wedgePhiVSkinPhi", sd1.detectorIndex);
	histoman->getHisto2D(hname)->Fill(pd1.phi, sd1.phi);


	hname = Form("hSABRE%d_EDif",sd1.detectorIndex);
	histoman->getHisto1D(hname)->Fill(sd1.ringEnergy-sd1.wedgeEnergy);
	

	hname = Form("hSABRE%d_PixelEDif",sd1.detectorIndex);
	histoman->getHisto3D(hname)->Fill(sd1.wedge, sd1.ring, sd1.ringEnergy - sd1.wedgeEnergy);
	

	hname = Form("hSABRE%d_ESummaryPixels",sd1.detectorIndex);
	histoman->getHisto3D(hname)->Fill(sd1.wedge, sd1.ring, sd1.ringEnergy);

	hname = Form("hSABRE%d_ChannelHits",sd1.detectorIndex);
	histoman->getHisto1D(hname)->Fill(globalring);
	histoman->getHisto1D(hname)->Fill(globalwedge);

	hname = Form("hSABRE%d_ERingSummary",sd1.detectorIndex);
	histoman->getHisto1D(hname)->Fill(sd1.ringEnergy);


	hname = Form("hSABRE%d_EWedgeSummary",sd1.detectorIndex);
	histoman->getHisto1D(hname)->Fill(sd1.wedgeEnergy);


	hname = Form("hSABRE%d_ChannelESummary",sd1.detectorIndex);
	histoman->getHisto2D(hname)->Fill(globalring,sd1.ringEnergy);
	histoman->getHisto2D(hname)->Fill(globalwedge,sd1.wedgeEnergy);


	if((sd1.detectorIndex == 3) && (sd1.ring == 7 || sd1.ring == 8 || sd1.ring == 9)){
		hname = Form("hSABRE%d_Ring%dSummary",sd1.detectorIndex,sd1.ring);
		histoman->getHisto1D(hname)->Fill(sd1.ringEnergy);
	}

	//all sabre histograms:
	histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalring);
	histoman->getHisto1D("hSABRE_ChannelHits")->Fill(globalwedge);
	histoman->getHisto1D("hSABRE_RingChannelHits")->Fill(globalring);
	histoman->getHisto1D("hSABRE_WedgeChannelHits")->Fill(globalwedge);
	histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalring,sd1.ringEnergy);
	histoman->getHisto2D("hSABRE_ChannelESummary")->Fill(globalwedge,sd1.wedgeEnergy);

}

void plot3mc::FillStraggleHistos(double oldTheta, double oldPhi, double newTheta, double newPhi, double dTheta, double dPhi){
	histoman->getHisto2D("hNewTheta_vs_OldTheta")->Fill(oldTheta, newTheta);
	histoman->getHisto2D("hNewPhi_vs_OldPhi")->Fill(oldPhi, newPhi);
	histoman->getHisto1D("hdTheta")->Fill(dTheta);
	histoman->getHisto1D("hdPhi")->Fill(dPhi);
}

bool plot3mc::FillTH1D(const TString& histoname, double value){
	if(!histoman->getHisto1D(histoname)) return false;
	histoman->getHisto1D(histoname)->Fill(value);
	return true;
}

bool plot3mc::FillTH2D(const TString& histoname, double valuex, double valuey){
	if(!histoman->getHisto2D(histoname)) return false;
	histoman->getHisto2D(histoname)->Fill(valuex, valuey);
	return true;
}

bool plot3mc::FillTH3D(const TString& histoname, double valuex, double valuey, double valuez){
	if(!histoman->getHisto3D(histoname)) return false;
	histoman->getHisto3D(histoname)->Fill(valuex, valuey, valuez);
	return true;
}

bool plot3mc::ParsePhysData(const std::string& line, PHYSDATA& pd1, PHYSDATA& pd2, PHYSDATA& pd3, PHYSDATA& pd4){
	std::istringstream iss(line);
	(iss >> pd1.e >> pd1.theta >> pd1.phi >> pd2.e >> pd2.theta >> pd2.phi >> pd3.e >> pd3.theta >> pd3.phi >> pd4.e >> pd4.theta >> pd4.phi);
	if(pd1.phi < 0) pd1.phi += 360.;
	if(pd2.phi < 0) pd2.phi += 360.;
	if(pd3.phi < 0) pd3.phi += 360.;
	if(pd4.phi < 0) pd4.phi += 360.;
	return true;
}

bool plot3mc::ParseSABREData(const std::string& line, SABREDATA& sd){
	std::istringstream iss(line);
	int index, ring, wedge;
	double ringe, wedgee, x, y;
	iss >> index >> ring >> wedge >> ringe >> wedgee >> x >> y;
	// (iss >> sd.detectorIndex >> sd.ring >> sd.wedge >> sd.ringEnergy >> sd.wedgeEnergy >> sd.localx >> sd.localy);
	sd.detectorIndex = index%100;
	sd.particleIndex = index-sd.detectorIndex;
	sd.ringEnergy = ringe;
	sd.wedgeEnergy = wedgee;
	sd.localx = x;
	sd.localy = y;
	sd.ring = ring;
	sd.wedge = wedge;
	sd.theta = -666.;//update this after the fact with angle map
	sd.phi = -666.;//update this after the fact with angle map
	return true;
}

void plot3mc::ProcessTXTOutput(std::vector<std::string> txtoutput){

	std::vector<std::string> eventLines;
	for(const auto& line : txtoutput){
		if(line == "-1"){

			//1 = ejectile
			//2 = recoil (not measured since breaks up)
			//3 = breakup1
			//4 = breakup2
			PHYSDATA pd1, pd2, pd3, pd4;
			SABREDATA sd1, sd2, sd3, sd4;

			if(eventLines.size() == 1){
				//only kinematics line (no particles in SABRE)
				ParsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				FillKinematicsHistos(pd1,pd2,pd3,pd4);


			} else if(eventLines.size() == 2){
				//kinematics line plus one particle in SABRE
				ParsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				ParseSABREData(eventLines[1],sd1);
				sd1.theta = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].first;//wow this is ugly but it works
				sd1.phi = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].second;//likewise^

				FillKinematicsHistos(pd1,pd2,pd3,pd4);

				if(sd1.particleIndex == 100){
					FillSABREHistos(sd1, pd1);
				} else if(sd1.particleIndex == 200){
					FillSABREHistos(sd1, pd2);
				} else if(sd1.particleIndex == 300){
					FillSABREHistos(sd1, pd3);
				} else if(sd1.particleIndex == 400){
					FillSABREHistos(sd1, pd4);
				}

			} else if(eventLines.size() == 3){
				//kinematics line plus two particles in SABRE
				ParsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				ParseSABREData(eventLines[1],sd1);
				sd1.theta = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].first;//wow this is ugly but it works
				sd1.phi = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].second;//likewise^
				ParseSABREData(eventLines[2],sd2);
				sd2.theta = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first, sd2.wedge+offsets[sd2.detectorIndex].second}].first;//wow this is ugly but it works
				sd2.phi = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first, sd2.wedge+offsets[sd2.detectorIndex].second}].second;//likewise^

				FillKinematicsHistos(pd1,pd2,pd3,pd4);

				if(sd1.particleIndex == 100){
					FillSABREHistos(sd1, pd1);
				} else if(sd1.particleIndex == 200){
					FillSABREHistos(sd1, pd2);
				} else if(sd1.particleIndex == 300){
					FillSABREHistos(sd1, pd3);
				} else if(sd1.particleIndex == 400){
					FillSABREHistos(sd1, pd4);
				}

				if(sd2.particleIndex == 100){
					FillSABREHistos(sd2, pd1);
				} else if(sd2.particleIndex == 200){
					FillSABREHistos(sd2, pd2);
				} else if(sd2.particleIndex == 300){
					FillSABREHistos(sd2, pd3);
				} else if(sd2.particleIndex == 400){
					FillSABREHistos(sd2, pd4);
				}

			} else if(eventLines.size() == 4){
				//kinematics line plus three particles in SABRE
				ParsePhysData(eventLines[0],pd1,pd2,pd3,pd4);
				ParseSABREData(eventLines[1],sd1);
				sd1.theta = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].first;//wow this is ugly but it works
				sd1.phi = sabre_thetaphimap[{sd1.ring+offsets[sd1.detectorIndex].first, sd1.wedge+offsets[sd1.detectorIndex].second}].second;//likewise^
				ParseSABREData(eventLines[2],sd2);
				sd2.theta = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first, sd2.wedge+offsets[sd2.detectorIndex].second}].first;//wow this is ugly but it works
				sd2.phi = sabre_thetaphimap[{sd2.ring+offsets[sd2.detectorIndex].first, sd2.wedge+offsets[sd2.detectorIndex].second}].second;//likewise^
				ParseSABREData(eventLines[3],sd3);
				sd3.theta = sabre_thetaphimap[{sd3.ring+offsets[sd3.detectorIndex].first, sd3.wedge+offsets[sd3.detectorIndex].second}].first;//wow this is ugly but it works
				sd3.phi = sabre_thetaphimap[{sd3.ring+offsets[sd3.detectorIndex].first, sd3.wedge+offsets[sd3.detectorIndex].second}].second;//likewise^

				if(sd1.particleIndex == 100){
					FillSABREHistos(sd1, pd1);
				} else if(sd1.particleIndex == 200){
					FillSABREHistos(sd1, pd2);
				} else if(sd1.particleIndex == 300){
					FillSABREHistos(sd1, pd3);
				} else if(sd1.particleIndex == 400){
					FillSABREHistos(sd1, pd4);
				}

				if(sd2.particleIndex == 100){
					FillSABREHistos(sd2, pd1);
				} else if(sd2.particleIndex == 200){
					FillSABREHistos(sd2, pd2);
				} else if(sd2.particleIndex == 300){
					FillSABREHistos(sd2, pd3);
				} else if(sd2.particleIndex == 400){
					FillSABREHistos(sd2, pd4);
				}

				if(sd3.particleIndex == 100){
					FillSABREHistos(sd3, pd1);
				} else if(sd3.particleIndex == 200){
					FillSABREHistos(sd3, pd2);
				} else if(sd3.particleIndex == 300){
					FillSABREHistos(sd3, pd3);
				} else if(sd3.particleIndex == 400){
					FillSABREHistos(sd3, pd4);
				}

			} else {
				std::cerr << "Warning: Unexpected eventLines.size() = " << eventLines.size() << "\n";
			}

			eventLines.clear();

		} else {

			eventLines.push_back(line);

		}
	}

}

void plot3mc::ProcessTXTOutput(const std::string& outputLines){

	std::vector<std::string> txtoutput = splitLines(outputLines);

	ProcessTXTOutput(txtoutput);
}

void plot3mc::SaveAndWrite(){
	UpdateHistoAxes();
	histoman->WriteAll(true);
}

double plot3mc::getSPSEnergy(double kinEnMeV, double sigmaMeV){
	double smearedE;

	smearedE = rand_gen->Gaus(kinEnMeV,sigmaMeV);

	return smearedE;
}

double plot3mc::calculateSPS_ExE(double spsE, double spsTheta, double spsPhi){//we should rethink this...hard coded for 6Li but also only used to fill explicitly 6Li ExE histos...
	TLorentzVector beam, target, ejectile, recoil;
	double smearedSPSE = getSPSEnergy(spsE);
	double beamEnergy = 7.5;
	beam.SetPxPyPzE(0.,0.,sqrt(2*masstable->GetMassMeV("He",3)*beamEnergy),beamEnergy+masstable->GetMassMeV("He",3));
	target.SetPxPyPzE(0.,0.,0.,masstable->GetMassMeV("Li",7));
	double pej = sqrt(2*smearedSPSE*masstable->GetMassMeV("He",4));
	ejectile.SetPxPyPzE(pej*sin(DEGRAD*spsTheta)*cos(DEGRAD*spsPhi), pej*sin(DEGRAD*spsTheta)*sin(DEGRAD*spsPhi), pej*cos(DEGRAD*spsTheta),smearedSPSE+masstable->GetMassMeV("He",4));
	recoil = beam + target - ejectile;

	return (recoil.M() - masstable->GetMassMeV("Li",6));
}

void plot3mc::readSingleAngleMap(std::ifstream& infile, std::map<std::pair<int,int>,std::pair<double,double>>& map){
	std::string header;
	getline(infile,header);

	int ringChannel, wedgeChannel;
	double theta, phi;
	while(infile >> ringChannel >> wedgeChannel >> theta >> phi){
		map[{ringChannel,wedgeChannel}] = {theta,phi};
	}
}

std::map<std::pair<int,int>, std::pair<double,double>> plot3mc::readAngleMaps(){
	const std::vector<std::string> filenames = {
		"anglemaps/SABRE0_phi306_anglemap.txt",
		"anglemaps/SABRE1_phi18_anglemap.txt",
		"anglemaps/SABRE2_phi234_anglemap.txt",
		"anglemaps/SABRE3_phi162_anglemap.txt",
		"anglemaps/SABRE4_phi90_anglemap.txt"
	};

	std::map<std::pair<int,int>, std::pair<double,double>> retmap;

	for(const auto& filename : filenames){
		std::ifstream infile(filename);
		if(!infile.is_open()){
			std::cerr << "Error: Failed to open file " << filename << "\n";
			continue;
		}
		readSingleAngleMap(infile,retmap);
		infile.close();
	}

	return retmap;
}

void plot3mc::UpdateHistoAxes(){

		for(int i=0; i<5; i++){
		TString hname = Form("hSABRE%d_hitsMapLocal",i);
		histoman->getHisto2D(hname)->GetXaxis()->SetTitle("+x <--------------------------------------------------> -x");
		histoman->getHisto2D(hname)->GetYaxis()->SetTitle("+y <--------------------------------------------------> -y");
		histoman->getHisto2D(hname)->GetXaxis()->CenterTitle();
		histoman->getHisto2D(hname)->GetYaxis()->CenterTitle();
	}

		std::vector<TString> histonames = {
		"hSABRE0_ringEVSkinE",
		"hSABRE1_ringEVSkinE",
		"hSABRE2_ringEVSkinE",
		"hSABRE3_ringEVSkinE",
		"hSABRE4_ringEVSkinE",
		"hSABRE0_wedgeEVSkinE",
		"hSABRE1_wedgeEVSkinE",
		"hSABRE2_wedgeEVSkinE",
		"hSABRE3_wedgeEVSkinE",
		"hSABRE4_wedgeEVSkinE",
		"hSABRE0_ringThetaVSkinTheta",
		"hSABRE1_ringThetaVSkinTheta",
		"hSABRE2_ringThetaVSkinTheta",
		"hSABRE3_ringThetaVSkinTheta",
		"hSABRE4_ringThetaVSkinTheta",
		"hSABRE0_wedgePhiVSkinPhi",
		"hSABRE1_wedgePhiVSkinPhi",
		"hSABRE2_wedgePhiVSkinPhi",
		"hSABRE3_wedgePhiVSkinPhi",
		"hSABRE4_wedgePhiVSkinPhi",
	};

	for(const auto& name : histonames){
		histoman->getHisto2D(name)->GetXaxis()->SetTitle("SABRESim");
		histoman->getHisto2D(name)->GetYaxis()->SetTitle("Kin2mc");
		histoman->getHisto2D(name)->GetXaxis()->CenterTitle();
		histoman->getHisto2D(name)->GetYaxis()->CenterTitle();

		//histoman->getHisto2D(name)->SetOption("SQUARE");
	}

	histonames = {
		"h_SABRETheta_vs_kinTheta",
		"h_SABREPhi_vs_kinPhi",
		"h_SABRERingE_vs_kinE",
		"h_SABREWedgeE_vs_kinE"
	};

	for(const auto& name : histonames){
		histoman->getHisto2D(name)->GetXaxis()->SetTitle("Kin2mc");
		histoman->getHisto2D(name)->GetYaxis()->SetTitle("SABREsim");
		histoman->getHisto2D(name)->GetXaxis()->CenterTitle();
		histoman->getHisto2D(name)->GetYaxis()->CenterTitle();
	}

	histonames = {
		"hSABRE0_EDif",
		"hSABRE1_EDif",
		"hSABRE2_EDif",
		"hSABRE3_EDif",
		"hSABRE4_EDif"
	};

	for(const auto& name : histonames){
		histoman->getHisto1D(name)->GetXaxis()->SetTitle("E_{ring} - E_{wedge}, any pixel");
		histoman->getHisto1D(name)->GetXaxis()->CenterTitle();
	}

	histonames = {
		"hSABRE0_PixelEDif",
		"hSABRE1_PixelEDif",
		"hSABRE2_PixelEDif",
		"hSABRE3_PixelEDif",
		"hSABRE4_PixelEDif"
	};

	for(const auto& name : histonames){
		histoman->getHisto3D(name)->GetXaxis()->SetTitle("Ring Index");
		histoman->getHisto3D(name)->GetYaxis()->SetTitle("Wedge Index");
		histoman->getHisto3D(name)->GetZaxis()->SetTitle("E_{ring} - E_{wedge} for pixel (ring,wedge)");
		histoman->getHisto3D(name)->GetXaxis()->CenterTitle();
		histoman->getHisto3D(name)->GetYaxis()->CenterTitle();
		histoman->getHisto3D(name)->GetZaxis()->CenterTitle();
	}

	histonames = {
		"hSABRE0_ESummaryPixels",
		"hSABRE1_ESummaryPixels",
		"hSABRE2_ESummaryPixels",
		"hSABRE3_ESummaryPixels",
		"hSABRE4_ESummaryPixels"
	};

	for(const auto& name : histonames){
		histoman->getHisto3D(name)->GetXaxis()->SetTitle("Wedge Index");
		histoman->getHisto3D(name)->GetYaxis()->SetTitle("Ring Index");
		histoman->getHisto3D(name)->GetZaxis()->SetTitle("ESummary By Pixel (using ERing)");
		histoman->getHisto3D(name)->GetXaxis()->CenterTitle();
		histoman->getHisto3D(name)->GetYaxis()->CenterTitle();
		histoman->getHisto3D(name)->GetZaxis()->CenterTitle();
	}

}
#include "invmass_mult3.h"
#include <iostream>
#include <vector>
#include <array>


InvMass_Mult3::InvMass_Mult3()
	: outfile(nullptr), intermediateMass(0), recoilMass(0), intermediateEx(0), intermediateExGate(0){

		permNames = {"012", "021", "102", "120", "201", "210", "allCases"};

		pMap = {
			{"012",{0,1,2}},
			{"021",{0,2,1}},
			{"102",{1,0,2}},
			{"120",{1,2,0}},
			{"201",{2,0,1}},
			{"210",{2,1,0}}
		};
	}

InvMass_Mult3::~InvMass_Mult3(){
	if(outfile && !outfile->IsZombie()){
		outfile->Close();
		delete outfile;
	}
}

void InvMass_Mult3::Init(const char* output_filename){
	
	outfile = new TFile(output_filename, "RECREATE");
	outtree = new TTree("InvMass_Mult3", "InvMass_Mult3");

	TString leaflist = "imIM/D:imEx/D:reconEx/D:imVCM/D:imKECM/D:imTHCM/D:imPHCM/D:"
					   "f1VCM/D:f1KECM/D:f1THCM/D:f1PHCM/D:"
					   "f2VCM/D:f2KECM/D:f2THCM/D:f2PHCM/D:"
					   "f3VCM/D:f3KECM/D:f3THCM/D:f3PHCM/D:"
					   "ecm1/D:ecm2/D";

	for(int i=0; i<6; i++){
		//create a branch for each permutation
		outtree->Branch(permNames[i], &caseResults[i], leaflist);
	}

	for(auto &cn : permNames){
		TDirectory *dir = outfile->mkdir(cn);
		dir->cd();

		// 1. Invariant Mass & Excitation Energy Histograms
		hMap[cn]["intermediateIM"]       = new TH1D(cn + "_intermediateIM", "Intermediate Invariant Mass;MeV/c^{2}", 29000, 4600, 7500);
		hMap[cn]["intermediateEx"]      = new TH1D(cn + "_intermediateEx", "Intermediate Ex;MeV", 525, -1, 20);
		hMap[cn]["ReconEx"]         = new TH1D(cn + "_ReconEx", "Recon Ex;MeV", 525, -1, 20);

		// 2. Kinematics Histograms (CM Frame)

		// Intermediate CM
		hMap[cn]["intermediatevcm"]      	   = new TH1D(cn + "_intermediatevcm", "Intermediate Velocity CM;c", 250, 0, 0.25);
		hMap[cn]["intermediatekecm"]     	   = new TH1D(cn + "_intermediatekecm", "Intermediate KE CM;MeV", 500, 0, 5);
		hMap[cn]["intermediatethetacm"]  	   = new TH1D(cn + "_intermediatethetacm", "Intermediate Theta CM;deg", 36, 0, 180);
		hMap[cn]["intermediatephicm"]    	   = new TH1D(cn + "_intermediatephicm", "Intermediate Phi CM;deg", 72, 0, 360);
		hMap[cn]["intermediatethetacmvsphicm"]  = new TH2D(cn + "_intermediatethetacmvsphicm", "Intermediate ThetaCM vs PhiCM;deg;deg", 72, 0, 360, 36, 0, 180);
		hMap[cn]["intermediatevcm_delta"]             = new TH1D(cn + "_intermediatevcm_delta", "intermediate Velocity CM (meas - expect);c", 500, -0.25, 0.25);
		hMap[cn]["intermediatekecm_delta"]            = new TH1D(cn + "_intermediatekecm_delta", "intermediate KE CM (meas - expect);MeV", 1000, -5, 5);

		// frag1
		hMap[cn]["frag1vcm"]             	   = new TH1D(cn + "_frag1vcm", "frag1 Velocity CM;c", 250, 0, 0.25);
		hMap[cn]["frag1kecm"]            	   = new TH1D(cn + "_frag1kecm", "frag1 KE CM;MeV", 500, 0, 5);
		hMap[cn]["frag1thetacm"]        	   = new TH1D(cn + "_frag1thetacm", "frag1 Theta CM;deg", 36, 0, 180);
		hMap[cn]["frag1phicm"]           	   = new TH1D(cn + "_frag1phicm", "frag1 Phi CM;deg", 72, 0, 360);
		hMap[cn]["frag1thetacmvsphicm"]        = new TH2D(cn + "_frag1thetacmvsphicm", "frag1 ThetaCM vs PhiCM;deg;deg", 72, 0, 360, 36, 0, 180);
		hMap[cn]["frag1vcm_delta"]             = new TH1D(cn + "_frag1vcm_delta", "frag1 Velocity CM (meas - expect);c", 500, -0.25, 0.25);
		hMap[cn]["frag1kecm_delta"]            = new TH1D(cn + "_frag1kecm_delta", "frag1 KE CM;MeV (meas - expect)", 1000, -5, 5);

		// decay1
		hMap[cn]["ecm1"]            	   	   = new TH1D(cn + "_ecm1", "E_{cm} Decay 1;MeV", 600, -1, 5);
		hMap[cn]["ecm1_delta"]             	   = new TH1D(cn + "_ecm1_delta", "E_{cm} Decay 1 (meas - expect);MeV", 1000, -5, 5);
		hMap[cn]["intermediatevcmVSfrag1vcm"]  = new TH2D(cn + "_intermediatevcmVSfrag1vcm", "Intermediate Vcm VS frag1 Vcm", 250, 0, 0.25, 250, 0, 0.25);
		hMap[cn]["intermediatekecmVSfrag1kecm"]= new TH2D(cn + "_intermediatekecmVSfrag1kecm", "Intermediate KEcm VS frag1 KEcm", 500, 0, 5, 500, 0, 5);

		// frag2
		hMap[cn]["frag2vcm"]        	       = new TH1D(cn + "_frag2vcm", "frag2 Velocity CM;c", 250, 0, 0.25);
		hMap[cn]["frag2kecm"]        	       = new TH1D(cn + "_frag2kecm", "frag2 KE CM;MeV", 500, 0, 5);
		hMap[cn]["frag2thetacm"]      		   = new TH1D(cn + "_frag2thetacm", "frag2 Theta CM;deg", 36, 0, 180);
		hMap[cn]["frag2phicm"]        		   = new TH1D(cn + "_frag2phicm", "frag2 Phi CM;deg", 72, 0, 360);
		hMap[cn]["frag2thetacmvsphicm"]        = new TH2D(cn + "_frag2thetacmvsphicm", "frag2 ThetaCM vs PhiCM;deg;deg", 72, 0, 360, 36, 0, 180);
		hMap[cn]["frag2vcm_delta"]             = new TH1D(cn + "_frag2vcm_delta", "frag2 Velocity CM (meas - expect);c", 500, -0.25, 0.25);
		hMap[cn]["frag2kecm_delta"]            = new TH1D(cn + "_frag2kecm_delta", "frag2 KE CM;MeV (meas - expect)", 1000, -5, 5);

		// frag3
		hMap[cn]["frag3vcm"]           		   = new TH1D(cn + "_frag3vcm", "frag3 Velocity CM;c", 250, 0, 0.25);
		hMap[cn]["frag3kecm"]         	   	   = new TH1D(cn + "_frag3kecm", "frag3 KE CM;MeV", 500, 0, 5);
		hMap[cn]["frag3thetacm"]      	  	   = new TH1D(cn + "_frag3thetacm", "frag3 Theta CM;deg", 36, 0, 180);
		hMap[cn]["frag3phicm"]      	       = new TH1D(cn + "_frag3phicm", "frag3 Phi CM;deg", 72, 0, 360);
		hMap[cn]["frag3thetacmvsphicm"]        = new TH2D(cn + "_frag3thetacmvsphicm", "frag3 ThetaCM vs PhiCM;deg;deg", 72, 0, 360, 36, 0, 180);
		hMap[cn]["frag3vcm_delta"]             = new TH1D(cn + "_frag3vcm_delta", "frag3 Velocity CM (meas - expect);c", 500, -0.25, 0.25);
		hMap[cn]["frag3kecm_delta"]            = new TH1D(cn + "_frag3kecm_delta", "frag3 KE CM;MeV (meas - expect)", 1000, -5, 5);

		// decay2
		hMap[cn]["ecm2"]             	   	   = new TH1D(cn + "_ecm2", "E_{cm} Decay 2;MeV", 600, -1, 5);
		hMap[cn]["ecm2_delta"]             	   = new TH1D(cn + "_ecm2_delta", "E_{cm} Decay 2 (meas - expect);MeV", 1000, -5, 5);
		hMap[cn]["frag2vcmVSfrag3vcm"]  	   = new TH2D(cn + "_frag2vcmVSfrag3vcm", "frag2 Vcm VS frag3 Vcm", 250, 0, 0.25, 250, 0, 0.25);
		hMap[cn]["frag2kecmVSfrag3kecm"]	   = new TH2D(cn + "_frag2kecmVSfrag3kecm", "frag2 KEcm VS frag3 KEcm", 500, 0, 5, 500, 0, 5);

		// Sequential Decay Energies
		hMap[cn]["ecm1VSecm2"]				   = new TH2D(cn + "_ecm1VSecm2", "ECM1 vs ECM2;E_{CM} Decay 2 (MeV); E_{CM} Decay 1 (MeV)", 600, -1, 5, 600, -1, 5);
		hMap[cn]["ecm1deltaVSecm2delta"]	   = new TH2D(cn + "_ecm1deltaVSecm2delta", "ECM1 delta vs ECM2 delta;E_{cm} Decay 2 (meas - expect); E_{CM} Decay 1 (meas - expect)", 1000, -5, 5, 1000, -5, 5);


		//gated histograms (gated on Ex of intermediate/intermediate)
		hMap[cn]["ReconEx_gated"]   = new TH1D(cn + "_ReconEx_gated", "Gated Recon Ex;MeV", 525, -1, 20);//gated on DaughterEx!

		hMap[cn]["intermediatevcm_gated"]      = new TH1D(cn + "_intermediatevcm_gated", "Intermediate Velocity CM;c", 250, 0, 0.25);
		hMap[cn]["intermediatekecm_gated"]     = new TH1D(cn + "_intermediatekecm_gated", "Intermediate KE CM;MeV", 500, 0, 5);
		hMap[cn]["intermediatethetacm_gated"]  = new TH1D(cn + "_intermediatethetacm_gated", "Intermediate Theta CM;deg", 36, 0, 180);
		hMap[cn]["intermediatephicm_gated"]    = new TH1D(cn + "_intermediatephicm_gated", "Intermediate Phi CM;deg", 72, 0, 360);
		hMap[cn]["intermediatethetacmvsphicm_gated"]  = new TH2D(cn + "_intermediatethetacmvsphicm_gated", "intermediate ThetaCM vs PhiCM;deg;deg", 72, 0, 360, 36, 0, 180);
		hMap[cn]["intermediatevcm_delta_gated"]             = new TH1D(cn + "_intermediatevcm_delta_gated", "intermediate Velocity CM (meas - expect);c", 500, -0.25, 0.25);
		hMap[cn]["intermediatekecm_delta_gated"]            = new TH1D(cn + "_intermediatekecm_delta_gated", "intermediate KE CM;MeV (meas - expect)", 1000, -5, 5);

		hMap[cn]["frag1vcm_gated"]             = new TH1D(cn + "_frag1vcm_gated", "frag1 Velocity CM;c", 250, 0, 0.25);
		hMap[cn]["frag1kecm_gated"]            = new TH1D(cn + "_frag1kecm_gated", "frag1 KE CM;MeV", 500, 0, 5);
		hMap[cn]["frag1thetacm_gated"]         = new TH1D(cn + "_frag1thetacm_gated", "frag1 Theta CM;deg", 36, 0, 180);
		hMap[cn]["frag1phicm_gated"]           = new TH1D(cn + "_frag1phicm_gated", "frag1 Phi CM;deg", 72, 0, 360);
		hMap[cn]["frag1thetacmvsphicm_gated"]  = new TH2D(cn + "_frag1thetacmvsphicm_gated", "frag1 ThetaCM vs PhiCM;deg;deg", 72, 0, 360, 36, 0, 180);
		hMap[cn]["frag1vcm_delta_gated"]       = new TH1D(cn + "_frag1vcm_delta_gated", "frag1 Velocity CM (meas - expect);c", 500, -0.25, 0.25);
		hMap[cn]["frag1kecm_delta_gated"]      = new TH1D(cn + "_frag1kecm_delta_gated", "frag1 KE CM;MeV (meas - expect)", 1000, -5, 5);

		hMap[cn]["intermediatevcmVSfrag1vcm_gated"]  = new TH2D(cn + "_intermediatevcmVSfrag1vcm_gated", "Intermediate Vcm VS frag1 Vcm", 250, 0, 0.25, 250, 0, 0.25);
		hMap[cn]["intermediatekecmVSfrag1kecm_gated"]= new TH2D(cn + "_intermediatekecmVSfrag1kecm_gated", "Intermediate KEcm VS frag1 KEcm", 500, 0, 5, 500, 0, 5);

		hMap[cn]["frag2vcm_gated"]             = new TH1D(cn + "_frag2vcm_gated", "frag2 Velocity CM;c", 250, 0, 0.25);
		hMap[cn]["frag2kecm_gated"]            = new TH1D(cn + "_frag2kecm_gated", "frag2 KE CM;MeV", 500, 0, 5);
		hMap[cn]["frag2thetacm_gated"]         = new TH1D(cn + "_frag2thetacm_gated", "frag2 Theta CM;deg", 36, 0, 180);
		hMap[cn]["frag2phicm_gated"]           = new TH1D(cn + "_frag2phicm_gated", "frag2 Phi CM;deg", 72, 0, 360);
		hMap[cn]["frag2thetacmvsphicm_gated"]  = new TH2D(cn + "_frag2thetacmvsphicm_gated", "frag2 ThetaCM vs PhiCM;deg;deg", 72, 0, 360, 36, 0, 180);
		hMap[cn]["frag2vcm_delta_gated"]       = new TH1D(cn + "_frag2vcm_delta_gated", "frag2 Velocity CM (meas - expect);c", 500, -0.25, 0.25);
		hMap[cn]["frag2kecm_delta_gated"]      = new TH1D(cn + "_frag2kecm_delta_gated", "frag2 KE CM;MeV (meas - expect)", 1000, -5, 5);

		hMap[cn]["frag3vcm_gated"]             = new TH1D(cn + "_frag3vcm_gated", "frag3 Velocity CM;c", 250, 0, 0.25);
		hMap[cn]["frag3kecm_gated"]            = new TH1D(cn + "_frag3kecm_gated", "frag3 KE CM;MeV", 500, 0, 5);
		hMap[cn]["frag3thetacm_gated"]         = new TH1D(cn + "_frag3thetacm_gated", "frag3 Theta CM;deg", 36, 0, 180);
		hMap[cn]["frag3phicm_gated"]           = new TH1D(cn + "_frag3phicm_gated", "frag3 Phi CM;deg", 72, 0, 360);
		hMap[cn]["frag3thetacmvsphicm_gated"]  = new TH2D(cn + "_frag3thetacmvsphicm_gated", "frag3 ThetaCM vs PhiCM;deg;deg", 72, 0, 360, 36, 0, 180);
		hMap[cn]["frag3vcm_delta_gated"]       = new TH1D(cn + "_frag3vcm_delta_gated", "frag3 Velocity CM (meas - expect);c", 500, -0.25, 0.25);
		hMap[cn]["frag3kecm_delta_gated"]      = new TH1D(cn + "_frag3kecm_delta_gated", "frag3 KE CM;MeV (meas - expect)", 1000, -5, 5);

		hMap[cn]["frag2vcmVSfrag3vcm_gated"]  	   = new TH2D(cn + "_frag2vcmVSfrag3vcm_gated", "frag2 Vcm VS frag3 Vcm", 250, 0, 0.25, 250, 0, 0.25);
		hMap[cn]["frag2kecmVSfrag3kecm_gated"]	   = new TH2D(cn + "_frag2kecmVSfrag3kecm_gated", "frag2 KEcm VS frag3 KEcm", 500, 0, 5, 500, 0, 5);



		hMap[cn]["ecm1_gated"]             	   = new TH1D(cn + "_ecm1_gated", "E_{cm} Decay 1;MeV", 600, -1, 5);
		hMap[cn]["ecm2_gated"]             	   = new TH1D(cn + "_ecm2_gated", "E_{cm} Decay 2;MeV", 600, -1, 5);
		hMap[cn]["ecm1VSecm2_gated"]				   = new TH2D(cn + "_ecm1VSecm2_gated", "ECM1 vs ECM2;E_{CM} Decay 2 (MeV); E_{CM} Decay 1 (MeV)", 600, -1, 5, 600, -1, 5);
		hMap[cn]["ecm1_delta_gated"]           = new TH1D(cn + "_ecm1_delta_gated", "E_{cm} Decay 1 (meas - expect);MeV", 1000, -5, 5);
		hMap[cn]["ecm2_delta_gated"]           = new TH1D(cn + "_ecm2_delta_gated", "E_{cm} Decay 2 (meas - expect);MeV", 1000, -5, 5);
		hMap[cn]["ecm1deltaVSecm2delta_gated"]	   = new TH2D(cn + "_ecm1deltaVSecm2delta_gated", "ECM1 delta vs ECM2 delta;E_{cm} Decay 2 (meas - expect); E_{CM} Decay 1 (meas - expect)", 1000, -5, 5, 1000, -5, 5);
		


		outfile->cd();
	}
}

void InvMass_Mult3::SetHypothesis(const Hypothesis4& hypo){
	hypothesis = hypo;

	for(int i=0; i<3; i++) masses[i] = hypo.masses[i];

	recoilMass = hypo.mass_recoil;
	intermediateMass = hypo.mass_intermediate;
	intermediateEx = hypo.intermediateEx;
	intermediateExGate = hypo.intermediateExGate;
	recoilEx = hypo.recoilEx;

	SetExpectedCMValues();
}

//AnalyzeEvent assumes theta, phi in degrees and E in MeV
//update this to return reconstructed recoil excitation energy in MeV
std::array<double,6> InvMass_Mult3::AnalyzeEvent(double E[3], double theta[3], double phi[3]){

	ClearEventResults();

	std::array<double,6> recoilExs;

	int permIndex = 0;
	for(auto const& [name, p] : pMap){

		TLorentzVector lv[3];
		int indices[3] = {p.i, p.j, p.k};

		//construct 4vectors for the hits based on current permutation (given in indices!)
		for(int n=0; n<3; n++){
			int hitindex = indices[n];
			double mass = masses[n];
			double mom = std::sqrt(2*mass*E[hitindex]);

			lv[n].SetPxPyPzE(
					mom*std::sin(theta[hitindex]*DEGRAD)*std::cos(phi[hitindex]*DEGRAD),
					mom*std::sin(theta[hitindex]*DEGRAD)*std::sin(phi[hitindex]*DEGRAD),
					mom*std::cos(theta[hitindex]*DEGRAD),
					E[hitindex] + mass
				);

		}

		//reconstruct relevant particles:
		TLorentzVector intermediate = lv[1] + lv[2];
		TLorentzVector recoil = lv[0] + lv[1] + lv[2];

		TLorentzVector frag1 = lv[0];
		TLorentzVector frag2 = lv[1];
		TLorentzVector frag3 = lv[2];

		//calculate excitation energy:
		double intermediateEx = intermediate.M() - intermediateMass;
		double Ex = recoil.M() - recoilMass;
		recoilExs[permIndex] = Ex;

		caseResults[permIndex].intermediateIM = intermediate.M();
		caseResults[permIndex].intermediateEx = intermediateEx;
		caseResults[permIndex].reconEx = Ex;

		TVector3 boost1 = -recoil.BoostVector();
		TVector3 boost2 = -intermediate.BoostVector();

		//begin analysis of first decay step: recoil -> frag1 + intermediate
		//boost the lab-measured intermediate and frag1 into the frame of the recoil:
		intermediate.Boost(boost1);
		frag1.Boost(boost1);

		double intermediatevcm = ((1/intermediate.Energy())*intermediate.Vect()).Mag();
		double intermediatekecm = 0.5*intermediateMass*intermediatevcm*intermediatevcm;
		double intermediatethetacm = RADDEG*std::acos(intermediate.Vect().Z()/intermediate.Vect().Mag());
		double intermediatephicm = RADDEG*std::atan2(intermediate.Vect().Y(), intermediate.Vect().X());
		if(intermediatephicm < 0) intermediatephicm += 360.;

		caseResults[permIndex].intermediatevcm = intermediatevcm;
		caseResults[permIndex].intermediatekecm = intermediatekecm;
		caseResults[permIndex].intermediatethetacm = intermediatethetacm;
		caseResults[permIndex].intermediatephicm = intermediatephicm;

		double frag1vcm = ((1/frag1.Energy())*frag1.Vect()).Mag();
		double frag1kecm = 0.5*masses[0]*frag1vcm*frag1vcm;
		double frag1thetacm = RADDEG*std::acos(frag1.Vect().Z()/frag1.Vect().Mag());
		double frag1phicm = RADDEG*std::atan2(frag1.Vect().Y(), frag1.Vect().X());
		if(frag1phicm < 0) frag1phicm += 360.;

		caseResults[permIndex].frag1vcm = frag1vcm;
		caseResults[permIndex].frag1kecm = frag1kecm;
		caseResults[permIndex].frag1thetacm = frag1thetacm;
		caseResults[permIndex].frag1phicm = frag1phicm;

		//determine ecm1:
		double ecm1 = intermediatekecm + frag1kecm;
		caseResults[permIndex].ecm1 = ecm1;


		//begin analysis of second decay step: intermediate -> frag2 + frag3
		//boost the lab-measured frag2 and frag3 into the frame of the intermediate:
		frag2.Boost(boost2);
		frag3.Boost(boost2);

		double frag2vcm = ((1/frag2.Energy())*frag2.Vect()).Mag();
		double frag2kecm = 0.5*masses[1]*frag2vcm*frag2vcm;
		double frag2thetacm = RADDEG*std::acos(frag2.Vect().Z()/frag2.Vect().Mag());
		double frag2phicm = RADDEG*std::atan2(frag2.Vect().Y(), frag2.Vect().X());
		if(frag2phicm < 0) frag2phicm += 360.;

		caseResults[permIndex].frag2vcm = frag2vcm;
		caseResults[permIndex].frag2kecm = frag2kecm;
		caseResults[permIndex].frag2thetacm = frag2thetacm;
		caseResults[permIndex].frag2phicm = frag2phicm;

		double frag3vcm = ((1/frag3.Energy())*frag3.Vect()).Mag();
		double frag3kecm = 0.5*masses[2]*frag3vcm*frag3vcm;
		double frag3thetacm = RADDEG*std::acos(frag3.Vect().Z()/frag3.Vect().Mag());
		double frag3phicm = RADDEG*std::atan2(frag3.Vect().Y(), frag3.Vect().X());
		if(frag3phicm < 0) frag3phicm += 360.;

		caseResults[permIndex].frag3vcm = frag3vcm;
		caseResults[permIndex].frag3kecm = frag3kecm;
		caseResults[permIndex].frag3thetacm = frag3thetacm;
		caseResults[permIndex].frag3phicm = frag3phicm;

		//determine ecm2:
		double ecm2 = frag2kecm + frag3kecm;
		caseResults[permIndex].ecm2 = ecm2;

		//increment permIndex here!
		permIndex += 1;

	}

	return recoilExs;
}

//FillEventHistograms calls FillSelectCaseHistograms for 0-5. Good for filling histograms regardless of results.
void InvMass_Mult3::FillEventHistograms(){
	for(int i=0; i<6; i++) FillSelectCaseHistograms(i);

	if(outtree) outtree->Fill();
}

//FillSelectCaseHistograms fills a single case histograms for caseNum = 0-5. Good for checking results conditionals before filling histograms.
void InvMass_Mult3::FillSelectCaseHistograms(int caseNum){
	//caseNum should be 0-5
	if(caseNum < 0 || caseNum > 5){
		std::cerr << "caseNum out of range (caseNum = " << caseNum << ", shuold be 0-5)\n";
		return;
	}

	hMap[permNames.at(caseNum)]["intermediateIM"]->Fill(caseResults[caseNum].intermediateIM);
	hMap[permNames.at(caseNum)]["intermediateEx"]->Fill(caseResults[caseNum].intermediateEx);
	hMap[permNames.at(caseNum)]["ReconEx"]->Fill(caseResults[caseNum].reconEx);
	hMap["allCases"]["intermediateIM"]->Fill(caseResults[caseNum].intermediateIM);
	hMap["allCases"]["intermediateEx"]->Fill(caseResults[caseNum].intermediateEx);
	hMap["allCases"]["ReconEx"]->Fill(caseResults[caseNum].reconEx);

	hMap[permNames.at(caseNum)]["intermediatevcm"]->Fill(caseResults[caseNum].intermediatevcm);
	hMap[permNames.at(caseNum)]["intermediatekecm"]->Fill(caseResults[caseNum].intermediatekecm);
	hMap[permNames.at(caseNum)]["intermediatevcm_delta"]->Fill(caseResults[caseNum].intermediatevcm - expectedCMValues.vcm_intermediate);
	hMap[permNames.at(caseNum)]["intermediatekecm_delta"]->Fill(caseResults[caseNum].intermediatekecm - expectedCMValues.kecm_intermediate);
	hMap[permNames.at(caseNum)]["intermediatethetacm"]->Fill(caseResults[caseNum].intermediatethetacm);
	hMap[permNames.at(caseNum)]["intermediatephicm"]->Fill(caseResults[caseNum].intermediatephicm);
	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["intermediatethetacmvsphicm"])->Fill(caseResults[caseNum].intermediatephicm, caseResults[caseNum].intermediatethetacm);
	hMap["allCases"]["intermediatevcm"]->Fill(caseResults[caseNum].intermediatevcm);
	hMap["allCases"]["intermediatekecm"]->Fill(caseResults[caseNum].intermediatekecm);
	hMap["allCases"]["intermediatevcm_delta"]->Fill(caseResults[caseNum].intermediatevcm - expectedCMValues.vcm_intermediate);
	hMap["allCases"]["intermediatekecm_delta"]->Fill(caseResults[caseNum].intermediatekecm - expectedCMValues.kecm_intermediate);
	hMap["allCases"]["intermediatethetacm"]->Fill(caseResults[caseNum].intermediatethetacm);
	hMap["allCases"]["intermediatephicm"]->Fill(caseResults[caseNum].intermediatephicm);
	dynamic_cast<TH2D*>(hMap["allCases"]["intermediatethetacmvsphicm"])->Fill(caseResults[caseNum].intermediatephicm, caseResults[caseNum].intermediatethetacm);

	hMap[permNames.at(caseNum)]["frag1vcm"]->Fill(caseResults[caseNum].frag1vcm);
	hMap[permNames.at(caseNum)]["frag1kecm"]->Fill(caseResults[caseNum].frag1kecm);
	hMap[permNames.at(caseNum)]["frag1vcm_delta"]->Fill(caseResults[caseNum].frag1vcm - expectedCMValues.vcm_frag1);
	hMap[permNames.at(caseNum)]["frag1kecm_delta"]->Fill(caseResults[caseNum].frag1kecm - expectedCMValues.kecm_frag1);
	hMap[permNames.at(caseNum)]["frag1thetacm"]->Fill(caseResults[caseNum].frag1thetacm);
	hMap[permNames.at(caseNum)]["frag1phicm"]->Fill(caseResults[caseNum].frag1phicm);
	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag1thetacmvsphicm"])->Fill(caseResults[caseNum].frag1phicm, caseResults[caseNum].frag1thetacm);
	hMap["allCases"]["frag1vcm"]->Fill(caseResults[caseNum].frag1vcm);
	hMap["allCases"]["frag1kecm"]->Fill(caseResults[caseNum].frag1kecm);
	hMap["allCases"]["frag1vcm_delta"]->Fill(caseResults[caseNum].frag1vcm - expectedCMValues.vcm_frag1);
	hMap["allCases"]["frag1kecm_delta"]->Fill(caseResults[caseNum].frag1kecm - expectedCMValues.kecm_frag1);
	hMap["allCases"]["frag1thetacm"]->Fill(caseResults[caseNum].frag1thetacm);
	hMap["allCases"]["frag1phicm"]->Fill(caseResults[caseNum].frag1phicm);
	dynamic_cast<TH2D*>(hMap["allCases"]["frag1thetacmvsphicm"])->Fill(caseResults[caseNum].frag1phicm, caseResults[caseNum].frag1thetacm);

	hMap[permNames.at(caseNum)]["ecm1"]->Fill(caseResults[caseNum].ecm1);
	hMap[permNames.at(caseNum)]["ecm1_delta"]->Fill(caseResults[caseNum].ecm1 - expectedCMValues.Ecm1);
	hMap["allCases"]["ecm1"]->Fill(caseResults[caseNum].ecm1);
	hMap["allCases"]["ecm1_delta"]->Fill(caseResults[caseNum].ecm1 - expectedCMValues.Ecm1);

	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["intermediatevcmVSfrag1vcm"])->Fill(caseResults[caseNum].frag1vcm, caseResults[caseNum].intermediatevcm);
	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["intermediatekecmVSfrag1kecm"])->Fill(caseResults[caseNum].frag1kecm, caseResults[caseNum].intermediatekecm);
	dynamic_cast<TH2D*>(hMap["allCases"]["intermediatevcmVSfrag1vcm"])->Fill(caseResults[caseNum].frag1vcm, caseResults[caseNum].intermediatevcm);
	dynamic_cast<TH2D*>(hMap["allCases"]["intermediatekecmVSfrag1kecm"])->Fill(caseResults[caseNum].frag1kecm, caseResults[caseNum].intermediatekecm);

	hMap[permNames.at(caseNum)]["frag2vcm"]->Fill(caseResults[caseNum].frag2vcm);
	hMap[permNames.at(caseNum)]["frag2kecm"]->Fill(caseResults[caseNum].frag2kecm);
	hMap[permNames.at(caseNum)]["frag2vcm_delta"]->Fill(caseResults[caseNum].frag2vcm - expectedCMValues.vcm_frag2);
	hMap[permNames.at(caseNum)]["frag2kecm_delta"]->Fill(caseResults[caseNum].frag2kecm - expectedCMValues.kecm_frag2);
	hMap[permNames.at(caseNum)]["frag2thetacm"]->Fill(caseResults[caseNum].frag2thetacm);
	hMap[permNames.at(caseNum)]["frag2phicm"]->Fill(caseResults[caseNum].frag2phicm);
	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag2thetacmvsphicm"])->Fill(caseResults[caseNum].frag2phicm, caseResults[caseNum].frag2thetacm);
	hMap["allCases"]["frag2vcm"]->Fill(caseResults[caseNum].frag2vcm);
	hMap["allCases"]["frag2kecm"]->Fill(caseResults[caseNum].frag2kecm);
	hMap["allCases"]["frag2vcm_delta"]->Fill(caseResults[caseNum].frag2vcm - expectedCMValues.vcm_frag2);
	hMap["allCases"]["frag2kecm_delta"]->Fill(caseResults[caseNum].frag2kecm - expectedCMValues.kecm_frag2);
	hMap["allCases"]["frag2thetacm"]->Fill(caseResults[caseNum].frag2thetacm);
	hMap["allCases"]["frag2phicm"]->Fill(caseResults[caseNum].frag2phicm);
	dynamic_cast<TH2D*>(hMap["allCases"]["frag2thetacmvsphicm"])->Fill(caseResults[caseNum].frag2phicm, caseResults[caseNum].frag2thetacm);

	hMap[permNames.at(caseNum)]["frag3vcm"]->Fill(caseResults[caseNum].frag3vcm);
	hMap[permNames.at(caseNum)]["frag3kecm"]->Fill(caseResults[caseNum].frag3kecm);
	hMap[permNames.at(caseNum)]["frag3vcm_delta"]->Fill(caseResults[caseNum].frag3vcm - expectedCMValues.vcm_frag3);
	hMap[permNames.at(caseNum)]["frag3kecm_delta"]->Fill(caseResults[caseNum].frag3kecm - expectedCMValues.kecm_frag3);
	hMap[permNames.at(caseNum)]["frag3thetacm"]->Fill(caseResults[caseNum].frag3thetacm);
	hMap[permNames.at(caseNum)]["frag3phicm"]->Fill(caseResults[caseNum].frag3phicm);
	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag3thetacmvsphicm"])->Fill(caseResults[caseNum].frag3phicm, caseResults[caseNum].frag3thetacm);
	hMap["allCases"]["frag3vcm"]->Fill(caseResults[caseNum].frag3vcm);
	hMap["allCases"]["frag3kecm"]->Fill(caseResults[caseNum].frag3kecm);
	hMap["allCases"]["frag3vcm_delta"]->Fill(caseResults[caseNum].frag3vcm - expectedCMValues.vcm_frag3);
	hMap["allCases"]["frag3kecm_delta"]->Fill(caseResults[caseNum].frag3kecm - expectedCMValues.kecm_frag3);
	hMap["allCases"]["frag3thetacm"]->Fill(caseResults[caseNum].frag3thetacm);
	hMap["allCases"]["frag3phicm"]->Fill(caseResults[caseNum].frag3phicm);
	dynamic_cast<TH2D*>(hMap["allCases"]["frag3thetacmvsphicm"])->Fill(caseResults[caseNum].frag3phicm, caseResults[caseNum].frag3thetacm);

	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag2vcmVSfrag3vcm"])->Fill(caseResults[caseNum].frag3vcm, caseResults[caseNum].frag2vcm);
	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag2kecmVSfrag3kecm"])->Fill(caseResults[caseNum].frag3kecm, caseResults[caseNum].frag2kecm);
	dynamic_cast<TH2D*>(hMap["allCases"]["frag2vcmVSfrag3vcm"])->Fill(caseResults[caseNum].frag3vcm, caseResults[caseNum].frag2vcm);
	dynamic_cast<TH2D*>(hMap["allCases"]["frag2kecmVSfrag3kecm"])->Fill(caseResults[caseNum].frag3kecm, caseResults[caseNum].frag2kecm);

	hMap[permNames.at(caseNum)]["ecm2"]->Fill(caseResults[caseNum].ecm2);
	hMap[permNames.at(caseNum)]["ecm2_delta"]->Fill(caseResults[caseNum].ecm2 - expectedCMValues.Ecm2);
	hMap["allCases"]["ecm2"]->Fill(caseResults[caseNum].ecm2);
	hMap["allCases"]["ecm2_delta"]->Fill(caseResults[caseNum].ecm2 - expectedCMValues.Ecm2);

	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["ecm1VSecm2"])->Fill(caseResults[caseNum].ecm2, caseResults[caseNum].ecm1);
	dynamic_cast<TH2D*>(hMap["allCases"]["ecm1VSecm2"])->Fill(caseResults[caseNum].ecm2, caseResults[caseNum].ecm1);

	dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["ecm1deltaVSecm2delta"])->Fill(caseResults[caseNum].ecm2 - expectedCMValues.Ecm2, caseResults[caseNum].ecm1 - expectedCMValues.Ecm1);
	dynamic_cast<TH2D*>(hMap["allCases"]["ecm1deltaVSecm2delta"])->Fill(caseResults[caseNum].ecm2 - expectedCMValues.Ecm2, caseResults[caseNum].ecm1 - expectedCMValues.Ecm1);

	if( std::abs(caseResults[caseNum].intermediateEx - intermediateEx) <= intermediateExGate ){
		hMap[permNames.at(caseNum)]["ReconEx_gated"]->Fill(caseResults[caseNum].reconEx);
		hMap["allCases"]["ReconEx_gated"]->Fill(caseResults[caseNum].reconEx);

		hMap[permNames.at(caseNum)]["intermediatevcm_gated"]->Fill(caseResults[caseNum].intermediatevcm);
		hMap[permNames.at(caseNum)]["intermediatekecm_gated"]->Fill(caseResults[caseNum].intermediatekecm);
		hMap[permNames.at(caseNum)]["intermediatevcm_delta_gated"]->Fill(caseResults[caseNum].intermediatevcm - expectedCMValues.vcm_intermediate);
		hMap[permNames.at(caseNum)]["intermediatekecm_delta_gated"]->Fill(caseResults[caseNum].intermediatekecm - expectedCMValues.kecm_intermediate);
		hMap[permNames.at(caseNum)]["intermediatethetacm_gated"]->Fill(caseResults[caseNum].intermediatethetacm);
		hMap[permNames.at(caseNum)]["intermediatephicm_gated"]->Fill(caseResults[caseNum].intermediatephicm);
		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["intermediatethetacmvsphicm_gated"])->Fill(caseResults[caseNum].intermediatephicm, caseResults[caseNum].intermediatethetacm);
		hMap["allCases"]["intermediatevcm_gated"]->Fill(caseResults[caseNum].intermediatevcm);
		hMap["allCases"]["intermediatekecm_gated"]->Fill(caseResults[caseNum].intermediatekecm);
		hMap["allCases"]["intermediatevcm_delta_gated"]->Fill(caseResults[caseNum].intermediatevcm - expectedCMValues.vcm_intermediate);
		hMap["allCases"]["intermediatekecm_delta_gated"]->Fill(caseResults[caseNum].intermediatekecm - expectedCMValues.kecm_intermediate);
		hMap["allCases"]["intermediatethetacm_gated"]->Fill(caseResults[caseNum].intermediatethetacm);
		hMap["allCases"]["intermediatephicm_gated"]->Fill(caseResults[caseNum].intermediatephicm);
		dynamic_cast<TH2D*>(hMap["allCases"]["intermediatethetacmvsphicm_gated"])->Fill(caseResults[caseNum].intermediatephicm, caseResults[caseNum].intermediatethetacm);

		hMap[permNames.at(caseNum)]["frag1vcm_gated"]->Fill(caseResults[caseNum].frag1vcm);
		hMap[permNames.at(caseNum)]["frag1kecm_gated"]->Fill(caseResults[caseNum].frag1kecm);
		hMap[permNames.at(caseNum)]["frag1vcm_delta_gated"]->Fill(caseResults[caseNum].frag1vcm - expectedCMValues.vcm_frag1);
		hMap[permNames.at(caseNum)]["frag1kecm_delta_gated"]->Fill(caseResults[caseNum].frag1kecm - expectedCMValues.kecm_frag1);
		hMap[permNames.at(caseNum)]["frag1thetacm_gated"]->Fill(caseResults[caseNum].frag1thetacm);
		hMap[permNames.at(caseNum)]["frag1phicm_gated"]->Fill(caseResults[caseNum].frag1phicm);
		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag1thetacmvsphicm_gated"])->Fill(caseResults[caseNum].frag1phicm, caseResults[caseNum].frag1thetacm);
		hMap["allCases"]["frag1vcm_gated"]->Fill(caseResults[caseNum].frag1vcm);
		hMap["allCases"]["frag1kecm_gated"]->Fill(caseResults[caseNum].frag1kecm);
		hMap["allCases"]["frag1vcm_delta_gated"]->Fill(caseResults[caseNum].frag1vcm - expectedCMValues.vcm_frag1);
		hMap["allCases"]["frag1kecm_delta_gated"]->Fill(caseResults[caseNum].frag1kecm - expectedCMValues.kecm_frag1);
		hMap["allCases"]["frag1thetacm_gated"]->Fill(caseResults[caseNum].frag1thetacm);
		hMap["allCases"]["frag1phicm_gated"]->Fill(caseResults[caseNum].frag1phicm);
		dynamic_cast<TH2D*>(hMap["allCases"]["frag1thetacmvsphicm_gated"])->Fill(caseResults[caseNum].frag1phicm, caseResults[caseNum].frag1thetacm);

		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["intermediatevcmVSfrag1vcm_gated"])->Fill(caseResults[caseNum].frag1vcm, caseResults[caseNum].intermediatevcm);
		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["intermediatekecmVSfrag1kecm_gated"])->Fill(caseResults[caseNum].frag1kecm, caseResults[caseNum].intermediatekecm);
		dynamic_cast<TH2D*>(hMap["allCases"]["intermediatevcmVSfrag1vcm_gated"])->Fill(caseResults[caseNum].frag1vcm, caseResults[caseNum].intermediatevcm);
		dynamic_cast<TH2D*>(hMap["allCases"]["intermediatekecmVSfrag1kecm_gated"])->Fill(caseResults[caseNum].frag1kecm, caseResults[caseNum].intermediatekecm);

		hMap[permNames.at(caseNum)]["frag2vcm_gated"]->Fill(caseResults[caseNum].frag2vcm);
		hMap[permNames.at(caseNum)]["frag2kecm_gated"]->Fill(caseResults[caseNum].frag2kecm);
		hMap[permNames.at(caseNum)]["frag2vcm_delta_gated"]->Fill(caseResults[caseNum].frag2vcm - expectedCMValues.vcm_frag2);
		hMap[permNames.at(caseNum)]["frag2kecm_delta_gated"]->Fill(caseResults[caseNum].frag2kecm - expectedCMValues.kecm_frag2);
		hMap[permNames.at(caseNum)]["frag2thetacm_gated"]->Fill(caseResults[caseNum].frag2thetacm);
		hMap[permNames.at(caseNum)]["frag2phicm_gated"]->Fill(caseResults[caseNum].frag2phicm);
		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag2thetacmvsphicm_gated"])->Fill(caseResults[caseNum].frag2phicm, caseResults[caseNum].frag2thetacm);
		hMap["allCases"]["frag2vcm_gated"]->Fill(caseResults[caseNum].frag2vcm);
		hMap["allCases"]["frag2kecm_gated"]->Fill(caseResults[caseNum].frag2kecm);
		hMap["allCases"]["frag2vcm_delta_gated"]->Fill(caseResults[caseNum].frag2vcm - expectedCMValues.vcm_frag2);
		hMap["allCases"]["frag2kecm_delta_gated"]->Fill(caseResults[caseNum].frag2kecm - expectedCMValues.kecm_frag2);
		hMap["allCases"]["frag2thetacm_gated"]->Fill(caseResults[caseNum].frag2thetacm);
		hMap["allCases"]["frag2phicm_gated"]->Fill(caseResults[caseNum].frag2phicm);
		dynamic_cast<TH2D*>(hMap["allCases"]["frag2thetacmvsphicm_gated"])->Fill(caseResults[caseNum].frag2phicm, caseResults[caseNum].frag2thetacm);

		hMap[permNames.at(caseNum)]["frag3vcm_gated"]->Fill(caseResults[caseNum].frag3vcm);
		hMap[permNames.at(caseNum)]["frag3kecm_gated"]->Fill(caseResults[caseNum].frag3kecm);
		hMap[permNames.at(caseNum)]["frag3vcm_delta_gated"]->Fill(caseResults[caseNum].frag3vcm - expectedCMValues.vcm_frag3);
		hMap[permNames.at(caseNum)]["frag3kecm_delta_gated"]->Fill(caseResults[caseNum].frag3kecm - expectedCMValues.kecm_frag3);
		hMap[permNames.at(caseNum)]["frag3thetacm_gated"]->Fill(caseResults[caseNum].frag3thetacm);
		hMap[permNames.at(caseNum)]["frag3phicm_gated"]->Fill(caseResults[caseNum].frag3phicm);
		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag3thetacmvsphicm_gated"])->Fill(caseResults[caseNum].frag3phicm, caseResults[caseNum].frag3thetacm);
		hMap["allCases"]["frag3vcm_gated"]->Fill(caseResults[caseNum].frag3vcm);
		hMap["allCases"]["frag3kecm_gated"]->Fill(caseResults[caseNum].frag3kecm);
		hMap["allCases"]["frag3vcm_delta_gated"]->Fill(caseResults[caseNum].frag3vcm - expectedCMValues.vcm_frag3);
		hMap["allCases"]["frag3kecm_delta_gated"]->Fill(caseResults[caseNum].frag3kecm - expectedCMValues.kecm_frag3);
		hMap["allCases"]["frag3thetacm_gated"]->Fill(caseResults[caseNum].frag3thetacm);
		hMap["allCases"]["frag3phicm_gated"]->Fill(caseResults[caseNum].frag3phicm);
		dynamic_cast<TH2D*>(hMap["allCases"]["frag3thetacmvsphicm_gated"])->Fill(caseResults[caseNum].frag3phicm, caseResults[caseNum].frag3thetacm);

		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag2vcmVSfrag3vcm_gated"])->Fill(caseResults[caseNum].frag3vcm, caseResults[caseNum].frag2vcm);
		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["frag2kecmVSfrag3kecm_gated"])->Fill(caseResults[caseNum].frag3kecm, caseResults[caseNum].frag2kecm);
		dynamic_cast<TH2D*>(hMap["allCases"]["frag2vcmVSfrag3vcm_gated"])->Fill(caseResults[caseNum].frag3vcm, caseResults[caseNum].frag2vcm);
		dynamic_cast<TH2D*>(hMap["allCases"]["frag2kecmVSfrag3kecm_gated"])->Fill(caseResults[caseNum].frag3kecm, caseResults[caseNum].frag2kecm);

		hMap[permNames.at(caseNum)]["ecm1_gated"]->Fill(caseResults[caseNum].ecm1);
		hMap[permNames.at(caseNum)]["ecm1_delta_gated"]->Fill(caseResults[caseNum].ecm1 - expectedCMValues.Ecm1);
		hMap[permNames.at(caseNum)]["ecm2_gated"]->Fill(caseResults[caseNum].ecm2);
		hMap[permNames.at(caseNum)]["ecm2_delta_gated"]->Fill(caseResults[caseNum].ecm2 - expectedCMValues.Ecm2);
		hMap["allCases"]["ecm1_gated"]->Fill(caseResults[caseNum].ecm1);
		hMap["allCases"]["ecm1_delta_gated"]->Fill(caseResults[caseNum].ecm1 - expectedCMValues.Ecm1);
		hMap["allCases"]["ecm2_gated"]->Fill(caseResults[caseNum].ecm2);
		hMap["allCases"]["ecm2_delta_gated"]->Fill(caseResults[caseNum].ecm2 - expectedCMValues.Ecm2);

		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["ecm1VSecm2_gated"])->Fill(caseResults[caseNum].ecm2, caseResults[caseNum].ecm1);
		dynamic_cast<TH2D*>(hMap["allCases"]["ecm1VSecm2_gated"])->Fill(caseResults[caseNum].ecm2, caseResults[caseNum].ecm1);
		dynamic_cast<TH2D*>(hMap[permNames.at(caseNum)]["ecm1deltaVSecm2delta_gated"])->Fill(caseResults[caseNum].ecm2 - expectedCMValues.Ecm2, caseResults[caseNum].ecm1 - expectedCMValues.Ecm1);
		dynamic_cast<TH2D*>(hMap["allCases"]["ecm1deltaVSecm2delta_gated"])->Fill(caseResults[caseNum].ecm2 - expectedCMValues.Ecm2, caseResults[caseNum].ecm1 - expectedCMValues.Ecm1);

	}
}

void InvMass_Mult3::CloseAndWrite(){
	if(outfile && outfile->IsOpen()){
		outfile->Write();
		outfile->Close();
	}
}

void InvMass_Mult3::SetExpectedCMValues(){
	double m_recoil = recoilMass + recoilEx;
	double m_inter = intermediateMass + intermediateEx;

	//decay 1 constants:
	expectedCMValues.Ecm1 = m_recoil - masses[0] - m_inter;
	if(expectedCMValues.Ecm1 > 0){
		expectedCMValues.kecm_frag1 = expectedCMValues.Ecm1 * (m_inter / (masses[0] + m_inter));
		expectedCMValues.kecm_intermediate = expectedCMValues.Ecm1 * (masses[0] / (masses[0] + m_inter));

		expectedCMValues.vcm_frag1 = std::sqrt(2.0 * expectedCMValues.kecm_frag1 / masses[0]);
		expectedCMValues.vcm_intermediate = std::sqrt(2.0 * expectedCMValues.kecm_intermediate / m_inter);
	}

	//decay 2 constants
	expectedCMValues.Ecm2 = m_inter - masses[1] - masses[2];
	if(expectedCMValues.Ecm2 > 0){
		expectedCMValues.kecm_frag2 = expectedCMValues.Ecm2 * (masses[2] / (masses[1] + masses[2]));
		expectedCMValues.kecm_frag3 = expectedCMValues.Ecm2 * (masses[1] / (masses[1] + masses[2]));

		expectedCMValues.vcm_frag2 = std::sqrt(2.0 * expectedCMValues.kecm_frag2 / masses[1]);
		expectedCMValues.vcm_frag3 = std::sqrt(2.0 * expectedCMValues.kecm_frag3 / masses[2]);
	}

	std::cout << "Masses:\n\tfrag1 = " << masses[0] << "\tfrag2 = " << masses[1] << "\tfrag3 = " << masses[2] << "\n";
	std::cout << "\tRecoil = " << recoilMass << " + " << recoilEx << " = " << m_recoil << "\n";
	std::cout << "\tIntermediate = " << intermediateMass << " + " << intermediateEx << " = " << m_inter << "\n";

	std::cout << "Decay 1 constants:" << std::endl;
	std::cout << "\tEcm1 = " << expectedCMValues.Ecm1 << std::endl;
	std::cout << "\tVcm  frag1 = " << expectedCMValues.vcm_frag1  << "\tKEcm frag1 = " << expectedCMValues.kecm_frag1 << std::endl;
	std::cout << "\tVcm  intermediate = " << expectedCMValues.vcm_intermediate  << "\tKEcm intermediate = " << expectedCMValues.kecm_intermediate << std::endl;

	std::cout << "\nDecay 2 constants:" << std::endl;
	std::cout << "\tEcm2 = " << expectedCMValues.Ecm2 << std::endl;
	std::cout << "\tVcm frag2 = " << expectedCMValues.vcm_frag2 << "\tKEcm frag2 = " << expectedCMValues.kecm_frag2 << std::endl;
	std::cout << "\tVcm frag3 = " << expectedCMValues.vcm_frag3 << "\tKEcm frag3 = " << expectedCMValues.kecm_frag3 << std::endl;
}

void InvMass_Mult3::ClearEventResults(){
	for(int i=0; i<6; i++){
		caseResults[i].Reset();
	}
}
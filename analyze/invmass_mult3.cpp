#include "invmass_mult3.h"
#include <iostream>
#include <vector>
#include <array>


InvMass_Mult3::InvMass_Mult3()
	: outfile(nullptr), intermediateMass(0), recoilMass(0){

		permNames = {"012", "021", "102", "120", "210", "201", "allCases"};

		pMap = {
			{"012",{0,1,2}},
			{"021",{0,2,1}},
			{"102",{1,0,2}},
			{"120",{1,2,0}},
			{"210",{2,1,0}},
			{"201",{2,0,1}}
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

	for(auto &cn : permNames){
		TDirectory *dir = outfile->mkdir(cn);
		dir->cd();

		// 1. Invariant Mass & Excitation Energy Histograms
		hMap[cn]["intermediateIM"]       = new TH1D(cn + "_intermediateIM", "Intermediate Invariant Mass;MeV/c^{2}", 29000, 4600, 7500);
		hMap[cn]["intermediateEx"]      = new TH1D(cn + "_intermediateEx", "Intermediate Ex;MeV", 525, -1, 20);
		hMap[cn]["ReconEx"]         = new TH1D(cn + "_ReconEx", "Recon Ex;MeV", 525, -1, 20);
		hMap[cn]["ReconEx_gated"]   = new TH1D(cn + "_ReconEx_gated", "Gated Recon Ex;MeV", 525, -1, 20);

		// 2. Kinematics Histograms (CM Frame)
		// Particle frag1
		hMap[cn]["frag1vcm"]             = new TH1D(cn + "_frag1vcm", "frag1 Velocity CM;c", 100, 0, 1);
		hMap[cn]["frag1kecm"]            = new TH1D(cn + "_frag1kecm", "frag1 KE CM;MeV", 500, 0, 5);
		hMap[cn]["frag1thetacm"]         = new TH1D(cn + "_frag1thetacm", "frag1 Theta CM;deg", 36, 0, 180);
		hMap[cn]["frag1phicm"]           = new TH1D(cn + "_frag1phicm", "frag1 Phi CM;deg", 72, 0, 360);

		// Particle frag2
		hMap[cn]["frag2vcm"]            = new TH1D(cn + "_frag2vcm", "frag2 Velocity CM;c", 100, 0, 1);
		hMap[cn]["frag2kecm"]           = new TH1D(cn + "_frag2kecm", "frag2 KE CM;MeV", 500, 0, 5);
		hMap[cn]["frag2thetacm"]        = new TH1D(cn + "_frag2thetacm", "frag2 Theta CM;deg", 36, 0, 180);
		hMap[cn]["frag2phicm"]          = new TH1D(cn + "_frag2phicm", "frag2 Phi CM;deg", 72, 0, 360);

		// Particle frag3
		hMap[cn]["frag3vcm"]            = new TH1D(cn + "_frag3vcm", "frag3 Velocity CM;c", 100, 0, 1);
		hMap[cn]["frag3kecm"]           = new TH1D(cn + "_frag3kecm", "frag3 KE CM;MeV", 500, 0, 5);
		hMap[cn]["frag3thetacm"]        = new TH1D(cn + "_frag3thetacm", "frag3 Theta CM;deg", 36, 0, 180);
		hMap[cn]["frag3phicm"]          = new TH1D(cn + "_frag3phicm", "frag3 Phi CM;deg", 72, 0, 360);

		// Intermediate specific CM
		hMap[cn]["intermediatevcm"]      = new TH1D(cn + "_intermediatevcm", "Intermediate Velocity CM;c", 100, 0, 1);
		hMap[cn]["intermediatekecm"]     = new TH1D(cn + "_intermediatekecm", "Intermediate KE CM;MeV", 500, 0, 5);
		hMap[cn]["intermediatethetacm"]  = new TH1D(cn + "_intermediatethetacm", "Intermediate Theta CM;deg", 36, 0, 180);
		hMap[cn]["intermediatephicm"]    = new TH1D(cn + "_intermediatephicm", "Intermediate Phi CM;deg", 72, 0, 360);

		// Sequential Decay Energies
		hMap[cn]["ecm1"]             = new TH1D(cn + "_ecm1", "E_{cm} Decay 1;MeV", 200, 0, 20);
		hMap[cn]["ecm2"]             = new TH1D(cn + "_ecm2", "E_{cm} Decay 2;MeV", 200, 0, 20);

		outfile->cd();
	}
}

//SetMasses assumes mass_frag1/2/3 are MeV/c^2
void InvMass_Mult3::SetMasses(double mass_frag1, double mass_frag2, double mass_frag3, double mass_recoil, double mass_intermediate){
	masses[0] = mass_frag1;
	masses[1] = mass_frag2;
	masses[2] = mass_frag3;
	recoilMass = mass_recoil;
	intermediateMass = mass_intermediate;
}

//ProcessEvent assumes theta, phi in degrees and E in MeV
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
		double recoilEx = recoil.M() - recoilMass;
		recoilExs[permIndex] = recoilEx;

		caseResults[permIndex].intermediateIM = intermediate.M();
		caseResults[permIndex].intermediateEx = intermediateEx;
		caseResults[permIndex].reconEx = recoilEx;

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
		// hMap[name]["ecm1"]->Fill(ecm1);
		// hMap["allCases"]["ecm1"]->Fill(ecm1);
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
	hMap[permNames.at(caseNum)]["intermediatethetacm"]->Fill(caseResults[caseNum].intermediatethetacm);
	hMap[permNames.at(caseNum)]["intermediatephicm"]->Fill(caseResults[caseNum].intermediatephicm);
	hMap["allCases"]["intermediatevcm"]->Fill(caseResults[caseNum].intermediatevcm);
	hMap["allCases"]["intermediatekecm"]->Fill(caseResults[caseNum].intermediatekecm);
	hMap["allCases"]["intermediatethetacm"]->Fill(caseResults[caseNum].intermediatethetacm);
	hMap["allCases"]["intermediatephicm"]->Fill(caseResults[caseNum].intermediatephicm);

	hMap[permNames.at(caseNum)]["frag1vcm"]->Fill(caseResults[caseNum].frag1vcm);
	hMap[permNames.at(caseNum)]["frag1kecm"]->Fill(caseResults[caseNum].frag1kecm);
	hMap[permNames.at(caseNum)]["frag1thetacm"]->Fill(caseResults[caseNum].frag1thetacm);
	hMap[permNames.at(caseNum)]["frag1phicm"]->Fill(caseResults[caseNum].frag1phicm);
	hMap["allCases"]["frag1vcm"]->Fill(caseResults[caseNum].frag1vcm);
	hMap["allCases"]["frag1kecm"]->Fill(caseResults[caseNum].frag1kecm);
	hMap["allCases"]["frag1thetacm"]->Fill(caseResults[caseNum].frag1thetacm);
	hMap["allCases"]["frag1phicm"]->Fill(caseResults[caseNum].frag1phicm);

	hMap[permNames.at(caseNum)]["ecm1"]->Fill(caseResults[caseNum].ecm1);
	hMap["allCases"]["ecm1"]->Fill(caseResults[caseNum].ecm1);

	hMap[permNames.at(caseNum)]["frag2vcm"]->Fill(caseResults[caseNum].frag2vcm);
	hMap[permNames.at(caseNum)]["frag2kecm"]->Fill(caseResults[caseNum].frag2kecm);
	hMap[permNames.at(caseNum)]["frag2thetacm"]->Fill(caseResults[caseNum].frag2thetacm);
	hMap[permNames.at(caseNum)]["frag2phicm"]->Fill(caseResults[caseNum].frag2phicm);
	hMap["allCases"]["frag2vcm"]->Fill(caseResults[caseNum].frag2vcm);
	hMap["allCases"]["frag2kecm"]->Fill(caseResults[caseNum].frag2kecm);
	hMap["allCases"]["frag2thetacm"]->Fill(caseResults[caseNum].frag2thetacm);
	hMap["allCases"]["frag2phicm"]->Fill(caseResults[caseNum].frag2phicm);

	hMap[permNames.at(caseNum)]["frag3vcm"]->Fill(caseResults[caseNum].frag3vcm);
	hMap[permNames.at(caseNum)]["frag3kecm"]->Fill(caseResults[caseNum].frag3kecm);
	hMap[permNames.at(caseNum)]["frag3thetacm"]->Fill(caseResults[caseNum].frag3thetacm);
	hMap[permNames.at(caseNum)]["frag3phicm"]->Fill(caseResults[caseNum].frag3phicm);
	hMap["allCases"]["frag3vcm"]->Fill(caseResults[caseNum].frag3vcm);
	hMap["allCases"]["frag3kecm"]->Fill(caseResults[caseNum].frag3kecm);
	hMap["allCases"]["frag3thetacm"]->Fill(caseResults[caseNum].frag3thetacm);
	hMap["allCases"]["frag3phicm"]->Fill(caseResults[caseNum].frag3phicm);

	hMap[permNames.at(caseNum)]["ecm2"]->Fill(caseResults[caseNum].ecm2);
	hMap["allCases"]["ecm2"]->Fill(caseResults[caseNum].ecm2);
}

void InvMass_Mult3::CloseAndWrite(){
	if(outfile && outfile->IsOpen()){
		outfile->Write();
		outfile->Close();
	}
}

// void InvMass_Mult3::SetExpectedCMValues(){
// 	//update expectedCMValues to have expected CM values
// 	expectedCMValues.Ecm1 = recoilMass - masses[0] - intermediateMass;
// 	expectedCMValues.Ecm2 = intermediateMass - masses[1] - masses[2];

// 	vcm_frag1 = sqrt(2*);
// }

void InvMass_Mult3::ClearEventResults(){
	for(int i=0; i<6; i++){
		caseResults[i].Reset();
	}
}
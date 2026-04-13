#include "invmass_mult3.h"
#include <iostream>


InvMass_Mult3::InvMass_Mult3()
	: outfile(nullptr), daughterMass(0), recoilMass(0){

		caseNames = {"012", "021", "102", "120", "210", "201", "allCases"};

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

	for(auto &cn : caseNames){
		TDirectory *dir = outfile->mkdir(cn);
		dir->cd();

		// 1. Invariant Mass & Excitation Energy Histograms
		hMap[cn]["daughterIM"]       = new TH1D(cn + "_daughterIM", "Daughter Invariant Mass;MeV/c^{2}", 29000, 4600, 7500);
		hMap[cn]["daughterExE"]      = new TH1D(cn + "_daughterExE", "Daughter Ex;MeV", 525, -1, 20);
		hMap[cn]["ReconExE"]         = new TH1D(cn + "_ReconExE", "9B Recon Ex;MeV", 525, -1, 20);
		hMap[cn]["ReconExE_gated"]   = new TH1D(cn + "_ReconExE_gated", "Gated 9B Recon Ex;MeV", 525, -1, 20);

		// 2. Kinematics Histograms (CM Frame)
		// Particle bu1
		hMap[cn]["bu1vcm"]             = new TH1D(cn + "_bu1vcm", "bu1 Velocity CM;c", 100, 0, 1);
		hMap[cn]["bu1kecm"]            = new TH1D(cn + "_bu1kecm", "bu1 KE CM;MeV", 100, 0, 1);
		hMap[cn]["bu1thetacm"]         = new TH1D(cn + "_bu1thetacm", "bu1 Theta CM;deg", 180, 0, 180);
		hMap[cn]["bu1phicm"]           = new TH1D(cn + "_bu1phicm", "bu1 Phi CM;deg", 360, 0, 360);

		// Particle bu2
		hMap[cn]["bu2vcm"]            = new TH1D(cn + "_bu2vcm", "bu2 Velocity CM;c", 100, 0, 1);
		hMap[cn]["bu2kecm"]           = new TH1D(cn + "_bu2kecm", "bu2 KE CM;MeV", 100, 0, 1);
		hMap[cn]["bu2thetacm"]        = new TH1D(cn + "_bu2thetacm", "bu2 Theta CM;deg", 180, 0, 180);
		hMap[cn]["bu2phicm"]          = new TH1D(cn + "_bu2phicm", "bu2 Phi CM;deg", 360, 0, 360);

		// Particle bu3
		hMap[cn]["bu3vcm"]            = new TH1D(cn + "_bu3vcm", "bu3 Velocity CM;c", 100, 0, 1);
		hMap[cn]["bu3kecm"]           = new TH1D(cn + "_bu3kecm", "bu3 KE CM;MeV", 100, 0, 1);
		hMap[cn]["bu3thetacm"]        = new TH1D(cn + "_bu3thetacm", "bu3 Theta CM;deg", 180, 0, 180);
		hMap[cn]["bu3phicm"]          = new TH1D(cn + "_bu3phicm", "bu3 Phi CM;deg", 360, 0, 360);

		// Daughter specific CM
		hMap[cn]["daughtervcm"]      = new TH1D(cn + "_daughtervcm", "Daughter Velocity CM;c", 100, 0, 1);
		hMap[cn]["daughterkecm"]     = new TH1D(cn + "_daughterkecm", "Daughter KE CM;MeV", 100, 0, 1);
		hMap[cn]["daughterthetacm"]  = new TH1D(cn + "_daughterthetacm", "Daughter Theta CM;deg", 180, 0, 180);
		hMap[cn]["daughterphicm"]    = new TH1D(cn + "_daughterphicm", "Daughter Phi CM;deg", 360, 0, 360);

		// Sequential Decay Energies
		hMap[cn]["ecm1"]             = new TH1D(cn + "_ecm1", "E_{cm} Decay 1;MeV", 200, 0, 20);
		hMap[cn]["ecm2"]             = new TH1D(cn + "_ecm2", "E_{cm} Decay 2;MeV", 200, 0, 20);

		outfile->cd();
	}

}

//SetMasses assumes mass_bu1/2/3 are MeV/c^2
void InvMass_Mult3::SetMasses(double mass_bu1, double mass_bu2, double mass_bu3, double mass_recoil, double mass_daughter){
	masses[0] = mass_bu1;
	masses[1] = mass_bu2;
	masses[2] = mass_bu3;
	recoilMass = mass_recoil;
	daughterMass = mass_daughter;
}

//ProcessEvent assumes theta, phi in degrees and E in MeV
void InvMass_Mult3::AnalyzeEvent(double E[3], double theta[3], double phi[3]){

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
		TLorentzVector daughter = lv[1] + lv[2];
		TLorentzVector recoil = lv[0] + lv[1] + lv[2];

		TLorentzVector bu1 = lv[0];
		TLorentzVector bu2 = lv[1];
		TLorentzVector bu3 = lv[2];

		//calculate excitation energy:
		double daughterEx = daughter.M() - daughterMass;
		double recoilEx = recoil.M() - recoilMass;

		//fill relevant histograms:
		hMap[name]["daughterIM"]->Fill(daughter.M());
		hMap[name]["daughterExE"]->Fill(daughterEx);
		hMap[name]["ReconExE"]->Fill(recoilEx);
		hMap["allCases"]["daughterIM"]->Fill(daughter.M());
		hMap["allCases"]["daughterExE"]->Fill(daughterEx);
		hMap["allCases"]["ReconExE"]->Fill(recoilEx);

		TVector3 boost1 = -recoil.BoostVector();
		TVector3 boost2 = -daughter.BoostVector();

		//begin analysis of first decay step: recoil -> bu1 + daughter
		//boost the lab-measured daughter and bu1 into the frame of the recoil:
		daughter.Boost(boost1);
		bu1.Boost(boost1);

		double daughtervcm = ((1/daughter.Energy())*daughter.Vect()).Mag();
		double daughterkecm = 0.5*daughterMass*daughtervcm*daughtervcm;
		double daughterthetacm = RADDEG*std::acos(daughter.Vect().Z()/daughter.Vect().Mag());
		double daughterphicm = RADDEG*std::atan2(daughter.Vect().Y(), daughter.Vect().X());
		if(daughterphicm < 0) daughterphicm += 360.;

		hMap[name]["daughtervcm"]->Fill(daughtervcm);
		hMap[name]["daughterkecm"]->Fill(daughterkecm);
		hMap[name]["daughterthetacm"]->Fill(daughterthetacm);
		hMap[name]["daughterphicm"]->Fill(daughterphicm);

		hMap["allCases"]["daughtervcm"]->Fill(daughtervcm);
		hMap["allCases"]["daughterkecm"]->Fill(daughterkecm);
		hMap["allCases"]["daughterthetacm"]->Fill(daughterthetacm);
		hMap["allCases"]["daughterphicm"]->Fill(daughterphicm);

		double bu1vcm = ((1/bu1.Energy())*bu1.Vect()).Mag();
		double bu1kecm = 0.5*masses[0]*bu1vcm*bu1vcm;
		double bu1thetacm = RADDEG*std::acos(bu1.Vect().Z()/bu1.Vect().Mag());
		double bu1phicm = RADDEG*std::atan2(bu1.Vect().Y(), bu1.Vect().X());
		if(bu1phicm < 0) bu1phicm += 360.;

		hMap[name]["bu1vcm"]->Fill(bu1vcm);
		hMap[name]["bu1kecm"]->Fill(bu1kecm);
		hMap[name]["bu1thetacm"]->Fill(bu1thetacm);
		hMap[name]["bu1phicm"]->Fill(bu1phicm);

		hMap["allCases"]["bu1vcm"]->Fill(bu1vcm);
		hMap["allCases"]["bu1kecm"]->Fill(bu1kecm);
		hMap["allCases"]["bu1thetacm"]->Fill(bu1thetacm);
		hMap["allCases"]["bu1phicm"]->Fill(bu1phicm);

		//determine ecm1:
		double ecm1 = daughterkecm + bu1kecm;
		hMap[name]["ecm1"]->Fill(ecm1);
		hMap["allCases"]["ecm1"]->Fill(ecm1);


		//begin analysis of second decay step: daughter -> bu2 + bu3
		//boost the lab-measured bu2 and bu3 into the frame of the daughter:
		bu2.Boost(boost2);
		bu3.Boost(boost2);

		double bu2vcm = ((1/bu2.Energy())*bu2.Vect()).Mag();
		double bu2kecm = 0.5*masses[1]*bu2vcm*bu2vcm;
		double bu2thetacm = RADDEG*std::acos(bu2.Vect().Z()/bu2.Vect().Mag());
		double bu2phicm = RADDEG*std::atan2(bu2.Vect().Y(), bu2.Vect().X());
		if(bu2phicm < 0) bu2phicm += 360.;

		hMap[name]["bu2vcm"]->Fill(bu2vcm);
		hMap[name]["bu2kecm"]->Fill(bu2kecm);
		hMap[name]["bu2thetacm"]->Fill(bu2thetacm);
		hMap[name]["bu2phicm"]->Fill(bu2phicm);

		hMap["allCases"]["bu2vcm"]->Fill(bu2vcm);
		hMap["allCases"]["bu2kecm"]->Fill(bu2kecm);
		hMap["allCases"]["bu2thetacm"]->Fill(bu2thetacm);
		hMap["allCases"]["bu2phicm"]->Fill(bu2phicm);

		double bu3vcm = ((1/bu3.Energy())*bu3.Vect()).Mag();
		double bu3kecm = 0.5*masses[2]*bu3vcm*bu3vcm;
		double bu3thetacm = RADDEG*std::acos(bu3.Vect().Z()/bu3.Vect().Mag());
		double bu3phicm = RADDEG*std::atan2(bu3.Vect().Y(), bu3.Vect().X());
		if(bu3phicm < 0) bu3phicm += 360.;

		hMap[name]["bu3vcm"]->Fill(bu3vcm);
		hMap[name]["bu3kecm"]->Fill(bu3kecm);
		hMap[name]["bu3thetacm"]->Fill(bu3thetacm);
		hMap[name]["bu3phicm"]->Fill(bu3phicm);

		hMap["allCases"]["bu3vcm"]->Fill(bu3vcm);
		hMap["allCases"]["bu3kecm"]->Fill(bu3kecm);
		hMap["allCases"]["bu3thetacm"]->Fill(bu3thetacm);
		hMap["allCases"]["bu3phicm"]->Fill(bu3phicm);


		//determine ecm2:
		double ecm2 = bu2kecm + bu3kecm;
		hMap[name]["ecm2"]->Fill(ecm2);
		hMap["allCases"]["ecm2"]->Fill(ecm2);

	}

}

void InvMass_Mult3::CloseAndWrite(){
	if(outfile && !outfile->IsZombie()){
		outfile->Write();
		outfile->Close();
	}

}


#include "invmass_mult3.h"
#include <iostream>
#include <vector>
#include <array>


InvMass_Mult3::InvMass_Mult3()
	: outfile(nullptr), intermediateMass(0), recoilMass(0), intermediateEx(0), intermediateEmin(0), intermediateEmax(0){

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
					   "ecm1/D:ecm2/D:"
					   "exp_ecm1/D:exp_ecm2/D:exp_imVCM/D:exp_imKECM/D:exp_imTHCM/D:exp_imPHCM/D:"
					   "exp_f1VCM/D:exp_f1KECM/D:exp_f2VCM/D:exp_f2KECM/D:"
					   "exp_f3VCM/D:exp_f3KECM/D";

	for(int i=0; i<6; i++){
		//create a branch for each permutation
		outtree->Branch(permNames[i], &caseResults[i], leaflist);
	}

	for(auto &cn : permNames){
		TDirectory *permdir = outfile->mkdir(cn);
		//dir->cd();

		TDirectory *ungateDir = permdir->mkdir("ungated");
		TDirectory *gateDir = permdir->mkdir("gated");

		groups_ungated[cn] = new permHisto(cn + "_ungated", ungateDir);
		groups_gated[cn] = new permHisto(cn + "_gated", gateDir);

		outfile->cd();
	}
}

void InvMass_Mult3::SetHypothesis(const Hypothesis4& hypo){
	hypothesis = hypo;

	for(int i=0; i<3; i++) masses[i] = hypo.masses[i];

	recoilMass = hypo.mass_recoil;
	intermediateMass = hypo.mass_intermediate;
	//intermediateEx = hypo.intermediateEx;
	//intermediateExGate = hypo.intermediateExGate;
	//recoilEx = hypo.recoilEx;

	//SetExpectedCMValues();
}

//AnalyzeEvent assumes theta, phi in degrees and E in MeV
//update this to return reconstructed recoil excitation energy in MeV
std::array<double,6> InvMass_Mult3::AnalyzeEvent(double E[3], double theta[3], double phi[3], bool updateIntermediateEx){

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

		//let's check if we should update CM variables by using the IM-calcualted intermediate excitation energy:
		if(updateIntermediateEx){
			SetIntermediateEx(intermediateEx);
			//test for dynamic4:
			//SetRecoilEx(Ex);
		}
		caseResults[permIndex].expected = expectedCMValues;

		caseResults[permIndex].intermediateIM = intermediate.M();
		caseResults[permIndex].intermediateEx = intermediateEx;
		caseResults[permIndex].reconEx = Ex;

		TVector3 boost1 = -recoil.BoostVector();
		TVector3 boost2 = -intermediate.BoostVector();

		caseResults[permIndex].boost1[0] = boost1.X();
		caseResults[permIndex].boost1[1] = boost1.Y();
		caseResults[permIndex].boost1[2] = boost1.Z();

		caseResults[permIndex].boost2[0] = boost2.X();
		caseResults[permIndex].boost2[1] = boost2.Y();
		caseResults[permIndex].boost2[2] = boost2.Z();

		//begin analysis of first decay step: recoil -> frag1 + intermediate
		//boost the lab-measured intermediate and frag1 into the frame of the recoil:
		intermediate.Boost(boost1);
		frag1.Boost(boost1);

		double intermediatevcm = ((1/intermediate.Energy())*intermediate.Vect()).Mag();
		//double intermediatevcm = intermediate.BoostVector().Mag();
		double intermediatekecm = 0.5*intermediateMass*intermediatevcm*intermediatevcm;
		double intermediatethetacm = RADDEG*std::acos(intermediate.Vect().Z()/intermediate.Vect().Mag());
		double intermediatephicm = RADDEG*std::atan2(intermediate.Vect().Y(), intermediate.Vect().X());
		if(intermediatephicm < 0) intermediatephicm += 360.;

		caseResults[permIndex].intermediatevcm = intermediatevcm;
		caseResults[permIndex].intermediatekecm = intermediatekecm;
		caseResults[permIndex].intermediatethetacm = intermediatethetacm;
		caseResults[permIndex].intermediatephicm = intermediatephicm;
		for(int i=0; i<3; i++) caseResults[permIndex].intermediateComp[i] = (-intermediate.BoostVector())(i);

		double frag1vcm = ((1/frag1.Energy())*frag1.Vect()).Mag();
		//double frag1vcm = frag1.BoostVector().Mag();
		double frag1kecm = 0.5*masses[0]*frag1vcm*frag1vcm;
		double frag1thetacm = RADDEG*std::acos(frag1.Vect().Z()/frag1.Vect().Mag());
		double frag1phicm = RADDEG*std::atan2(frag1.Vect().Y(), frag1.Vect().X());
		if(frag1phicm < 0) frag1phicm += 360.;

		caseResults[permIndex].frag1vcm = frag1vcm;
		caseResults[permIndex].frag1kecm = frag1kecm;
		caseResults[permIndex].frag1thetacm = frag1thetacm;
		caseResults[permIndex].frag1phicm = frag1phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag1Comp[i] = (-frag1.BoostVector())(i);

		//determine ecm1:
		double ecm1 = intermediatekecm + frag1kecm + intermediateEx;
		caseResults[permIndex].ecm1 = ecm1;


		//begin analysis of second decay step: intermediate -> frag2 + frag3
		//boost the lab-measured frag2 and frag3 into the frame of the intermediate:
		frag2.Boost(boost2);
		frag3.Boost(boost2);

		double frag2vcm = ((1/frag2.Energy())*frag2.Vect()).Mag();
		//double frag2vcm = frag2.BoostVector().Mag();
		double frag2kecm = 0.5*masses[1]*frag2vcm*frag2vcm;
		double frag2thetacm = RADDEG*std::acos(frag2.Vect().Z()/frag2.Vect().Mag());
		double frag2phicm = RADDEG*std::atan2(frag2.Vect().Y(), frag2.Vect().X());
		if(frag2phicm < 0) frag2phicm += 360.;

		caseResults[permIndex].frag2vcm = frag2vcm;
		caseResults[permIndex].frag2kecm = frag2kecm;
		caseResults[permIndex].frag2thetacm = frag2thetacm;
		caseResults[permIndex].frag2phicm = frag2phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag2Comp[i] = (-frag2.BoostVector())(i);

		double frag3vcm = ((1/frag3.Energy())*frag3.Vect()).Mag();
		//double frag3vcm = frag3.BoostVector().Mag();
		double frag3kecm = 0.5*masses[2]*frag3vcm*frag3vcm;
		double frag3thetacm = RADDEG*std::acos(frag3.Vect().Z()/frag3.Vect().Mag());
		double frag3phicm = RADDEG*std::atan2(frag3.Vect().Y(), frag3.Vect().X());
		if(frag3phicm < 0) frag3phicm += 360.;

		caseResults[permIndex].frag3vcm = frag3vcm;
		caseResults[permIndex].frag3kecm = frag3kecm;
		caseResults[permIndex].frag3thetacm = frag3thetacm;
		caseResults[permIndex].frag3phicm = frag3phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag3Comp[i] = (-frag3.BoostVector())(i);

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

void InvMass_Mult3::FillSelectCaseHistograms(int caseNum){
	if(caseNum < 0 || caseNum > 5){
		std::cerr << "caseNum out of range (" << caseNum << ")\n";
		return;
	}

	TString cn = permNames.at(caseNum);
	auto& res = caseResults[caseNum];

	//lambdas to fill both specific permutation and "allCases"
	auto fillAll = [&](TString key, double x){
		groups_ungated[cn]->Fill(key, x);
		groups_ungated["allCases"]->Fill(key, x);
	};

	auto fillAll2D = [&](TString key, double x, double y){
		groups_ungated[cn]->Fill(key, x, y);
		groups_ungated["allCases"]->Fill(key, x, y);
	};

	//invariant mass and excitation energy histograms:
	fillAll("intermediateIM", res.intermediateIM);
	fillAll("intermediateEx", res.intermediateEx);
	fillAll("RecoilEx", res.reconEx);

	//intermediate CM
	fillAll("intermediatevcm_meas", res.intermediatevcm);
	fillAll("intermediatevcm_expect", res.expected.vcm_intermediate);
	fillAll("intermediatevcm_delta", res.intermediatevcm - res.expected.vcm_intermediate);
	fillAll2D("intermediatevcm_TransverseVSLongitudinal", std::abs(res.intermediateComp[2]), std::sqrt(res.intermediateComp[0]*res.intermediateComp[0] + res.intermediateComp[1]*res.intermediateComp[1]));

	fillAll("intermediatekecm_meas", res.intermediatekecm);
	fillAll("intermediatekecm_expect", res.expected.kecm_intermediate);
	fillAll("intermediatekecm_delta", res.intermediatekecm - res.expected.kecm_intermediate);


	fillAll("intermediatethetacm", res.intermediatethetacm);
	fillAll("intermediatephicm", res.intermediatephicm);
	fillAll2D("intermediatethetacmvsphicm", res.intermediatephicm, res.intermediatethetacm);
	fillAll2D("intermediatevcmVSthetacm",res.intermediatethetacm, res.intermediatevcm);
	fillAll2D("intermediatekecmVSthetacm",res.intermediatethetacm, res.intermediatekecm);

	//fill frag1-3
	auto fillFrag = [&](int i, double vcm, double comps[3], double kecm, double theta, double phi, double expV, double expK){
		TString f = Form("frag%d",i);
		fillAll(f+"vcm_meas",vcm);
		fillAll(f+"vcm_expect",expV);
		fillAll(f+"vcm_delta",vcm-expV);
		fillAll2D(f+"vcm_TransverseVSLongitudinal", std::abs(comps[2]), std::sqrt(comps[0]*comps[0] + comps[1]*comps[1]));

		fillAll(f+"kecm_meas",kecm);
		fillAll(f+"kecm_expect",expK);
		fillAll(f+"kecm_delta",kecm-expK);

		fillAll(f+"thetacm",theta);
		fillAll(f+"phicm",phi);
		fillAll2D(f+"thetacmvsphicm",phi,theta);

		fillAll2D(f+"vcmVSthetacm",theta,vcm);
		fillAll2D(f+"kecmVSthetacm",theta,kecm);
	};
	fillFrag(1, res.frag1vcm, res.frag1Comp, res.frag1kecm, res.frag1thetacm, res.frag1phicm, res.expected.vcm_frag1, res.expected.kecm_frag1);
	fillFrag(2, res.frag2vcm, res.frag2Comp, res.frag2kecm, res.frag2thetacm, res.frag2phicm, res.expected.vcm_frag2, res.expected.kecm_frag2);
	fillFrag(3, res.frag3vcm, res.frag3Comp, res.frag3kecm, res.frag3thetacm, res.frag3phicm, res.expected.vcm_frag3, res.expected.kecm_frag3);

	//decays
	fillAll("ecm1_meas", res.ecm1);
	fillAll("ecm1_expect", res.expected.Ecm1);
	fillAll("ecm1_delta", res.ecm1 - res.expected.Ecm1);
	fillAll2D("ecm1measVSintermediatethetacm", res.intermediatethetacm, res.ecm1);
	fillAll2D("ecm1measVSfrag1thetacm", res.frag1thetacm, res.ecm1);
	fillAll2D("ecm1measVSfrag2thetacm", res.frag2thetacm, res.ecm1);
	fillAll2D("ecm1measVSfrag3thetacm", res.frag3thetacm, res.ecm1);
	fillAll("decay1_VCM", std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1] + res.boost1[2]*res.boost1[2]));
	fillAll2D("decay1_VCM_TransverseVSLongitudinal", std::abs(res.boost1[2]), std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1]));
	fillAll("decay1_thetaCMsum", res.intermediatethetacm + res.frag1thetacm);
	fillAll("decay1_phiCMdiff", std::abs(res.intermediatephicm - res.frag1phicm));
	fillAll("ecm2_meas", res.ecm2);
	fillAll("ecm2_expect", res.expected.Ecm2);
	fillAll("ecm2_delta", res.ecm2 - res.expected.Ecm2);
	fillAll2D("ecm2measVSintermediatethetacm", res.intermediatethetacm, res.ecm2);
	fillAll2D("ecm2measVSfrag1thetacm", res.frag1thetacm, res.ecm2);
	fillAll2D("ecm2measVSfrag2thetacm", res.frag2thetacm, res.ecm2);
	fillAll2D("ecm2measVSfrag3thetacm", res.frag3thetacm, res.ecm2);
	fillAll("decay2_VCM", std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1] + res.boost2[2]*res.boost2[2]));
	fillAll2D("decay2_VCM_TransverseVSLongitudinal", std::abs(res.boost2[2]), std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1]));
	fillAll("decay2_thetaCMsum", res.frag2thetacm + res.frag3thetacm);
	fillAll("decay2_phiCMdiff", std::abs(res.frag2phicm - res.frag3phicm));

	fillAll2D("intermediatevcmVSfrag1vcm", res.frag1vcm, res.intermediatevcm);
	fillAll2D("intermediatekecmVSfrag1kecm", res.frag1kecm, res.intermediatekecm);
	fillAll2D("frag2vcmVSfrag3vcm", res.frag3vcm, res.frag2vcm);
	fillAll2D("frag2kecmVSfrag3kecm", res.frag3kecm, res.frag2kecm);
	fillAll2D("ecm1VSecm2", res.ecm2, res.ecm1);
	fillAll2D("ecm1deltaVSecm2delta", res.ecm2 - res.expected.Ecm2, res.ecm1 - res.expected.Ecm1);

	bool intermediateExCheck = (res.intermediateEx >= intermediateEmin && res.intermediateEx <= intermediateEmax);
	bool intermediateVcmCheck = (res.intermediatevcm >= 0.0092 && res.intermediatevcm <= 0.0102);
	bool frag1VcmCheck = (res.frag1vcm >= 0.0735 && res.frag1vcm <= 0.0795);

	//if( std::abs(res.intermediateEx - intermediateEx) <= intermediateExGate ){
	//if(res.intermediateEx >= intermediateEmin && res.intermediateEx <= intermediateEmax){
	if(intermediateExCheck){
		if(intermediateVcmCheck){
			//lambdas to fill both specific permutation and "allCases"
			auto fillGated = [&](TString key, double x){
				groups_gated[cn]->Fill(key, x);
				groups_gated["allCases"]->Fill(key, x);
			};

			auto fillGated2D = [&](TString key, double x, double y){
				groups_gated[cn]->Fill(key, x, y);
				groups_gated["allCases"]->Fill(key, x, y);
			};

			//invariant mass and excitation energy histograms:
			fillGated("intermediateIM", res.intermediateIM);
			fillGated("intermediateEx", res.intermediateEx);
			fillGated("RecoilEx", res.reconEx);

			//intermediate CM
			fillGated("intermediatevcm_meas", res.intermediatevcm);
			fillGated("intermediatevcm_expect", res.expected.vcm_intermediate);
			fillGated("intermediatevcm_delta", res.intermediatevcm - res.expected.vcm_intermediate);
			fillGated2D("intermediatevcm_TransverseVSLongitudinal", std::abs(res.intermediateComp[2]), std::sqrt(res.intermediateComp[0]*res.intermediateComp[0] + res.intermediateComp[1]*res.intermediateComp[1]));

			fillGated("intermediatekecm_meas", res.intermediatekecm);
			fillGated("intermediatekecm_expect", res.expected.kecm_intermediate);
			fillGated("intermediatekecm_delta", res.intermediatekecm - res.expected.kecm_intermediate);

			fillGated("intermediatethetacm", res.intermediatethetacm);
			fillGated("intermediatephicm", res.intermediatephicm);
			fillGated2D("intermediatethetacmvsphicm", res.intermediatephicm, res.intermediatethetacm);

			fillGated2D("intermediatevcmVSthetacm",res.intermediatethetacm, res.intermediatevcm);
			fillGated2D("intermediatekecmVSthetacm",res.intermediatethetacm, res.intermediatekecm);

			//fill frag1-3
			auto fillFragGated = [&](int i, double vcm, double comps[3], double kecm, double theta, double phi, double expV, double expK){
				TString f = Form("frag%d",i);
				fillGated(f+"vcm_meas",vcm);
				fillGated(f+"vcm_expect",expV);
				fillGated(f+"vcm_delta",vcm-expV);
				fillGated2D(f+"vcm_TransverseVSLongitudinal", std::abs(comps[2]), std::sqrt(comps[0]*comps[0] + comps[1]*comps[1]));

				fillGated(f+"kecm_meas",kecm);
				fillGated(f+"kecm_expect",expK);
				fillGated(f+"kecm_delta",kecm-expK);

				fillGated(f+"thetacm",theta);
				fillGated(f+"phicm",phi);
				fillGated2D(f+"thetacmvsphicm",phi,theta);

				fillGated2D(f+"vcmVSthetacm",theta,vcm);
				fillGated2D(f+"kecmVSthetacm",theta,kecm);
			};
			fillFragGated(1, res.frag1vcm, res.frag1Comp, res.frag1kecm, res.frag1thetacm, res.frag1phicm, res.expected.vcm_frag1, res.expected.kecm_frag1);
			fillFragGated(2, res.frag2vcm, res.frag2Comp, res.frag2kecm, res.frag2thetacm, res.frag2phicm, res.expected.vcm_frag2, res.expected.kecm_frag2);
			fillFragGated(3, res.frag3vcm, res.frag3Comp, res.frag3kecm, res.frag3thetacm, res.frag3phicm, res.expected.vcm_frag3, res.expected.kecm_frag3);

			//decays
			fillGated("ecm1_meas", res.ecm1);
			fillGated("ecm1_expect", res.expected.Ecm1);
			fillGated("ecm1_delta", res.ecm1 - res.expected.Ecm1);
			fillGated2D("ecm1measVSintermediatethetacm", res.intermediatethetacm, res.ecm1);
			fillGated2D("ecm1measVSfrag1thetacm", res.frag1thetacm, res.ecm1);
			fillGated2D("ecm1measVSfrag2thetacm", res.frag2thetacm, res.ecm1);
			fillGated2D("ecm1measVSfrag3thetacm", res.frag3thetacm, res.ecm1);
			fillGated("decay1_VCM", std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1] + res.boost1[2]*res.boost1[2]));
			fillGated2D("decay1_VCM_TransverseVSLongitudinal", std::abs(res.boost1[2]), std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1]));
			fillGated("decay1_thetaCMsum", res.intermediatethetacm + res.frag1thetacm);
			fillGated("decay1_phiCMdiff", std::abs(res.intermediatephicm - res.frag1phicm));
			fillGated("ecm2_meas", res.ecm2);
			fillGated("ecm2_expect", res.expected.Ecm2);
			fillGated("ecm2_delta", res.ecm2 - res.expected.Ecm2);
			fillGated2D("ecm2measVSintermediatethetacm", res.intermediatethetacm, res.ecm2);
			fillGated2D("ecm2measVSfrag1thetacm", res.frag1thetacm, res.ecm2);
			fillGated2D("ecm2measVSfrag2thetacm", res.frag2thetacm, res.ecm2);
			fillGated2D("ecm2measVSfrag3thetacm", res.frag3thetacm, res.ecm2);
			fillGated("decay2_VCM", std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1] + res.boost2[2]*res.boost2[2]));
			fillGated2D("decay2_VCM_TransverseVSLongitudinal", std::abs(res.boost2[2]), std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1]));
			fillGated("decay2_thetaCMsum", res.frag2thetacm + res.frag3thetacm);
			fillGated("decay2_phiCMdiff", std::abs(res.frag2phicm - res.frag3phicm));

			fillGated2D("intermediatevcmVSfrag1vcm", res.frag1vcm, res.intermediatevcm);
			fillGated2D("intermediatekecmVSfrag1kecm", res.frag1kecm, res.intermediatekecm);
			fillGated2D("frag2vcmVSfrag3vcm", res.frag3vcm, res.frag2vcm);
			fillGated2D("frag2kecmVSfrag3kecm", res.frag3kecm, res.frag2kecm);
			fillGated2D("ecm1VSecm2", res.ecm2, res.ecm1);
			fillGated2D("ecm1deltaVSecm2delta", res.ecm2 - res.expected.Ecm2, res.ecm1 - res.expected.Ecm1);
		}
	}
}

void InvMass_Mult3::CloseAndWrite(){
	if(outfile && outfile->IsOpen()){
		outfile->Write();
		outfile->Close();
	}
}

void InvMass_Mult3::SetExpectedCMValues(bool verbose){
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

	if(verbose){
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

}

void InvMass_Mult3::ClearEventResults(){
	for(int i=0; i<6; i++){
		caseResults[i].Reset();
	}
}
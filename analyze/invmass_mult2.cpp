#include "invmass_mult2.h"
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

InvMass_Mult2::InvMass_Mult2() : outfile(nullptr), recoilMass(0), recoilEx(0){

	permNames = {"01", "10", "allCases"};

	pMap = {
		{"01", {0,1}},
		{"10", {1,0}}
	};
}

InvMass_Mult2::~InvMass_Mult2(){
	if(outfile && !outfile->IsZombie()){
		outfile->Close();
		delete outfile;
	}
}

void InvMass_Mult2::Init(const char* output_filename){

	outfile = new TFile(output_filename, "RECREATE");
	outtree = new TTree("InvMass_Mult2", "InvMass_Mult2");

	TString leaflist = "recoilIM/D:reconEx/D:"
					   "f1VCM/D:f1KECM/D:f1THCM/D:f1PHCM/D:f1Comp[3]/D:"
					   "f2VCM/D:f2KECM/D:f2THCM/D:f2PHCM/D:f2Comp[3]/D:"
					   "ecm/D:boost[3]/D:relLabAngle_frag1frag2/D:"
					   "exp_ecm/D:"
					   "exp_f1VCM/D:exp_f1KECM/D:"
					   "exp_f2VCM/D:exp_f2KECM/D:"
					   "permPasses/O";

	for(int i=0; i<2; i++){
		outtree->Branch("_" + permNames[i], &caseResults[i], leaflist);
	}

	for(auto &pn : permNames){
		TDirectory *permdir = outfile->mkdir(pn);

		TDirectory *ungateDir = permdir->mkdir("ungated");
		TDirectory *gateDir = permdir->mkdir("gated");

		groups_ungated[pn] = new permHisto_mult2(pn + "_ungated", ungateDir);
		groups_gated[pn] = new permHisto_mult2(pn + "_gated", gateDir);

		outfile->cd();
	}

	hPermCounter = new TH1D("hPermCounter", "hPermCounter", 3, -0.5, 2.5);
	hPermCounter_gated = new TH1D("hPermCounter_gated", "hPermCounter_gated", 3, -0.5, 2.5);
	hSortedPermutations = new TH1D("hSortedPermutations", "Permutation Picks", 2, -0.5, 1.5);
	hSortedIMRecEx = new TH1D("hSortedIMRecEx", "Rec Ex (IM)", 300, -5, 7);

	const char* labels[2] = {"01", "10"};
	for(int i=0; i<2; i++){
		hSortedPermutations->GetXaxis()->SetBinLabel(i + 1, labels[i]);
	}
}

void InvMass_Mult2::SetHypothesis(const Hypothesis3& hypo){
	hypothesis = hypo;

	for(int i=0; i<2; i++) masses[i] = hypo.masses[i];
	recoilMass = hypo.mass_recoil;
}

std::array<double, 2> InvMass_Mult2::AnalyzeEvent(double E[2], double theta[2], double phi[2]){

	ClearEventResults();
	std::array<double, 2> recoilExs;

	int permIndex = 0;
	for(auto const& [name, p] : pMap){

		TLorentzVector lv[2];
		int indices[2] = {p.i, p.j};

		for(int n=0; n<2; n++){
			int hitindex = indices[n];
			double mass = masses[n];
			double mom = std::sqrt(2 * mass * E[hitindex]);

			lv[n].SetPxPyPzE(
				mom * std::sin(theta[hitindex]*DEGRAD) * std::cos(phi[hitindex]*DEGRAD),
				mom * std::sin(theta[hitindex]*DEGRAD) * std::sin(phi[hitindex]*DEGRAD),
				mom * std::cos(theta[hitindex]*DEGRAD),
				E[hitindex] + mass
			);
		}

		TLorentzVector recoil = lv[0] + lv[1];
		TLorentzVector frag1 = lv[0];
		TLorentzVector frag2 = lv[1];

		caseResults[permIndex].relLabAngle_frag1frag2 = frag1.Vect().Angle(frag2.Vect())*RADDEG;

		double Ex = recoil.M() - recoilMass;
		recoilExs[permIndex] = Ex;

		caseResults[permIndex].permName = name;
		caseResults[permIndex].expected = expectedCMValues;

		caseResults[permIndex].recoilIM = recoil.M();
		caseResults[permIndex].reconEx = Ex;

		TVector3 boostVector = -recoil.BoostVector();
		caseResults[permIndex].boost[0] = boostVector.X();
		caseResults[permIndex].boost[1] = boostVector.Y();
		caseResults[permIndex].boost[2] = boostVector.Z();

		frag1.Boost(boostVector);
		frag2.Boost(boostVector);

		double frag1vcm = ((1/frag1.Energy())*frag1.Vect()).Mag();
		double frag1kecm = 0.5*masses[0]*frag1vcm*frag1vcm;
		double frag1thetacm = RADDEG*std::acos(frag1.Vect().Z() / frag1.Vect().Mag());
		double frag1phicm = RADDEG*std::atan2(frag1.Vect().Y(), frag1.Vect().X());
		if(frag1phicm < 0.) frag1phicm += 360.;

		caseResults[permIndex].frag1vcm = frag1vcm;
		caseResults[permIndex].frag1kecm = frag1kecm;
		caseResults[permIndex].frag1thetacm = frag1thetacm;
		caseResults[permIndex].frag1phicm = frag1phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag1Comp[i] = (-frag1.BoostVector())(i);

		double frag2vcm = ((1/frag2.Energy())*frag2.Vect()).Mag();
		double frag2kecm = 0.5*masses[1]*frag2vcm*frag2vcm;
		double frag2thetacm = RADDEG*std::acos(frag2.Vect().Z() / frag2.Vect().Mag());
		double frag2phicm = RADDEG*std::atan2(frag2.Vect().Y(), frag2.Vect().X());
		if(frag2phicm < 0.) frag2phicm += 360.;

		caseResults[permIndex].frag2vcm = frag2vcm;
		caseResults[permIndex].frag2kecm = frag2kecm;
		caseResults[permIndex].frag2thetacm = frag2thetacm;
		caseResults[permIndex].frag2phicm = frag2phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag2Comp[i] = (-frag2.BoostVector())(i);

		double ecm = frag1kecm + frag2kecm;
		caseResults[permIndex].ecm = ecm;

		caseResults[permIndex].exp_ecm = expectedCMValues.Ecm;
		caseResults[permIndex].exp_f1VCM = expectedCMValues.vcm_frag1;
		caseResults[permIndex].exp_f1KECM = expectedCMValues.kecm_frag1;
		caseResults[permIndex].exp_f2VCM = expectedCMValues.vcm_frag2;
		caseResults[permIndex].exp_f2KECM = expectedCMValues.kecm_frag2;

		permIndex += 1;

	}

	//gate evaluation logic:
	for(int i=0; i<2; i++){
		bool gate1 = false, gate2 = false;
		auto& res = caseResults[i];

		switch(gate1_index){
			case NOCHECK_M2:
				gate1 = true;
				break;
			case FRAG1VCMCHECK_M2:
				gate1 = CheckGate1(res.frag1vcm);
				break;
			case FRAG2VCMCHECK_M2:
				gate1 = CheckGate1(res.frag2vcm);
				break;
			default:
				gate1 = false;
				break;
		}

		switch(gate2_index){
			case NOCHECK_M2:
				gate2 = true;
				break;
			case FRAG1VCMCHECK_M2:
				gate2 = CheckGate2(res.frag1vcm);
				break;
			case FRAG2VCMCHECK_M2:
				gate2 = CheckGate2(res.frag2vcm);
				break;
			default:
				gate2 = false;
				break;
		}

		if(gate1 && gate2) res.permPasses = true;
	}

	return recoilExs;
}

void InvMass_Mult2::FillTree(){
	if(outtree) outtree->Fill();
}

void InvMass_Mult2::FillEventHistograms(double SPS_Ex){
	for(int i=0; i<2; i++) FillSelectCaseHistograms(i, SPS_Ex);
}

void InvMass_Mult2::FillSelectCaseHistograms(int caseNum, double SPS_Ex) {
	if (caseNum < 0 || caseNum > 1) {
		std::cerr << "caseNum out of range (" << caseNum << ")\n";
		return;
	}

	TString cn = permNames.at(caseNum);
	auto& res = caseResults[caseNum];

	auto fillAll = [&](TString key, double x) {
		groups_ungated[cn]->Fill(key, x);
		groups_ungated["allCases"]->Fill(key, x);
	};

	auto fillAll2D = [&](TString key, double x, double y) {
		groups_ungated[cn]->Fill(key, x, y);
		groups_ungated["allCases"]->Fill(key, x, y);
	};

	fillAll("RecoilIM", res.recoilIM);
	fillAll("RecoilEx", res.reconEx);
	fillAll2D("RecoilEx_IMvsSPS", SPS_Ex, res.reconEx);
	fillAll("RecoilExDif", SPS_Ex - res.reconEx);

	auto fillFrag = [&](int i, double vcm, double comps[3], double kecm, double theta, double phi, double expV, double expK) {
		TString f = Form("frag%d", i);
		fillAll(f + "vcm_meas", vcm);
		//fillAll(f + "vcm_expect", expV);
		//fillAll(f + "vcm_delta", vcm - expV);
		fillAll2D(f + "vcm_TransverseVSLongitudinal", std::abs(comps[2]), std::sqrt(comps[0]*comps[0] + comps[1]*comps[1]));

		fillAll(f + "kecm_meas", kecm);
		//fillAll(f + "kecm_expect", expK);
		//fillAll(f + "kecm_delta", kecm - expK);

		fillAll(f + "thetacm", theta);
		fillAll(f + "phicm", phi);
		fillAll2D(f + "thetacmvsphicm", phi, theta);
		fillAll2D(f + "vcmVSthetacm", theta, vcm);
		fillAll2D(f + "kecmVSthetacm", theta, kecm);
    };

	fillFrag(1, res.frag1vcm, res.frag1Comp, res.frag1kecm, res.frag1thetacm, res.frag1phicm, res.expected.vcm_frag1, res.expected.kecm_frag1);
	fillFrag(2, res.frag2vcm, res.frag2Comp, res.frag2kecm, res.frag2thetacm, res.frag2phicm, res.expected.vcm_frag2, res.expected.kecm_frag2);

	fillAll("ecm_meas", res.ecm);
	//fillAll("ecm_expect", res.expected.Ecm);
	//fillAll("ecm_delta", res.ecm - res.expected.Ecm);
	fillAll2D("ecmmeasVSfrag1thetacm", res.frag1thetacm, res.ecm);
	fillAll2D("ecmmeasVSfrag2thetacm", res.frag2thetacm, res.ecm);

	fillAll("decay_VCM", std::sqrt(res.boost[0]*res.boost[0] + res.boost[1]*res.boost[1] + res.boost[2]*res.boost[2]));
	fillAll2D("decay_VCM_TransverseVSLongitudinal", std::abs(res.boost[2]), std::sqrt(res.boost[0]*res.boost[0] + res.boost[1]*res.boost[1]));
	fillAll("decay_thetaCMsum", res.frag1thetacm + res.frag2thetacm);
	fillAll("decay_phiCMdiff", std::abs(res.frag1phicm - res.frag2phicm));
	fillAll("decay_relLabAngle", res.relLabAngle_frag1frag2);

	fillAll2D("frag1vcmVSfrag2vcm", res.frag2vcm, res.frag1vcm);
	fillAll2D("frag1kecmVSfrag2kecm", res.frag2kecm, res.frag1kecm);
}

void InvMass_Mult2::FillGatedEventHistograms(double SPS_Ex){
	for(int i=0; i<2; i++) FillSelectGatedCaseHistograms(i, SPS_Ex);
}

void InvMass_Mult2::FillSelectGatedCaseHistograms(int caseNum, double SPS_Ex) {
	if (caseNum < 0 || caseNum > 1) {
		std::cerr << "caseNum out of range (" << caseNum << ")\n";
		return;
	}

	TString cn = permNames.at(caseNum);
	auto& res = caseResults[caseNum];

	if (res.permPasses) {
		auto fillGated = [&](TString key, double x) {
			groups_gated[cn]->Fill(key, x);
			groups_gated["allCases"]->Fill(key, x);
		};

		auto fillGated2D = [&](TString key, double x, double y) {
			groups_gated[cn]->Fill(key, x, y);
			groups_gated["allCases"]->Fill(key, x, y);
		};

		fillGated("RecoilIM", res.recoilIM);
		fillGated("RecoilEx", res.reconEx);
		fillGated2D("RecoilEx_IMvsSPS", SPS_Ex, res.reconEx);
		fillGated("RecoilExDif", SPS_Ex - res.reconEx);

		auto fillFragGated = [&](int i, double vcm, double comps[3], double kecm, double theta, double phi, double expV, double expK) {
			TString f = Form("frag%d", i);
			fillGated(f + "vcm_meas", vcm);
			//fillGated(f + "vcm_expect", expV);
			//fillGated(f + "vcm_delta", vcm - expV);
			fillGated2D(f + "vcm_TransverseVSLongitudinal", std::abs(comps[2]), std::sqrt(comps[0]*comps[0] + comps[1]*comps[1]));

			fillGated(f + "kecm_meas", kecm);
			//fillGated(f + "kecm_expect", expK);
			//fillGated(f + "kecm_delta", kecm - expK);

			fillGated(f + "thetacm", theta);
			fillGated(f + "phicm", phi);
			fillGated2D(f + "thetacmvsphicm", phi, theta);
			fillGated2D(f + "vcmVSthetacm", theta, vcm);
			fillGated2D(f + "kecmVSthetacm", theta, kecm);
		};

		fillFragGated(1, res.frag1vcm, res.frag1Comp, res.frag1kecm, res.frag1thetacm, res.frag1phicm, res.expected.vcm_frag1, res.expected.kecm_frag1);
		fillFragGated(2, res.frag2vcm, res.frag2Comp, res.frag2kecm, res.frag2thetacm, res.frag2phicm, res.expected.vcm_frag2, res.expected.kecm_frag2);

		fillGated("ecm_meas", res.ecm);
		//fillGated("ecm_expect", res.expected.Ecm);
		//fillGated("ecm_delta", res.ecm - res.expected.Ecm);
		fillGated2D("ecmmeasVSfrag1thetacm", res.frag1thetacm, res.ecm);
		fillGated2D("ecmmeasVSfrag2thetacm", res.frag2thetacm, res.ecm);

		fillGated("decay_VCM", std::sqrt(res.boost[0]*res.boost[0] + res.boost[1]*res.boost[1] + res.boost[2]*res.boost[2]));
		fillGated2D("decay_VCM_TransverseVSLongitudinal", std::abs(res.boost[2]), std::sqrt(res.boost[0]*res.boost[0] + res.boost[1]*res.boost[1]));
		fillGated("decay_thetaCMsum", res.frag1thetacm + res.frag2thetacm);
		fillGated("decay_phiCMdiff", std::abs(res.frag1phicm - res.frag2phicm));
		fillGated("decay_relLabAngle", res.relLabAngle_frag1frag2);

		fillGated2D("frag1vcmVSfrag2vcm", res.frag2vcm, res.frag1vcm);
		fillGated2D("frag1kecmVSfrag2kecm", res.frag2kecm, res.frag1kecm);
	}
}

void InvMass_Mult2::FillPermCounter(bool gated){
	if(gated){
		hPermCounter_gated->Fill(CountPermPasses());
	} else {
		hPermCounter->Fill(CountPermPasses());
	}
}

void InvMass_Mult2::FillSortedHisto(double SPS_Ex){
	Results sorted_caseResults[2];
	for(int i=0; i<2; i++) sorted_caseResults[i] = caseResults[i];

	std::sort(std::begin(sorted_caseResults), std::end(sorted_caseResults), [SPS_Ex](const Results& a, const Results& b) {
		return std::abs(SPS_Ex - a.reconEx) < std::abs(SPS_Ex - b.reconEx);
	});

	hSortedIMRecEx->Fill(sorted_caseResults[0].reconEx);

	int permIndex = -1;
	TString perm = sorted_caseResults[0].permName;
	if(perm == "01") permIndex = 0;
	else if(perm == "10") permIndex = 1;

	hSortedPermutations->Fill(permIndex);
}

void InvMass_Mult2::CloseAndWrite(){
	outfile->cd();
	hPermCounter->Write();
	hPermCounter_gated->Write();
	hSortedPermutations->Write();
	hSortedIMRecEx->Write();
	if(outfile && outfile->IsOpen()){
		outfile->Write();
		outfile->Close();
	}
}

void InvMass_Mult2::SetExpectedCMValues(bool verbose){
	double m_recoil = recoilMass + recoilEx;

	expectedCMValues.Ecm = m_recoil - masses[0] - masses[1];
	if(expectedCMValues.Ecm > 0){
		expectedCMValues.kecm_frag1 = expectedCMValues.Ecm * (masses[1] / (masses[0] + masses[1]));
		expectedCMValues.kecm_frag2 = expectedCMValues.Ecm * (masses[0] / (masses[0] + masses[1]));

		expectedCMValues.vcm_frag1 = std::sqrt(2.0 * expectedCMValues.kecm_frag1 / masses[0]);
		expectedCMValues.vcm_frag2 = std::sqrt(2.0 * expectedCMValues.kecm_frag2 / masses[1]);
	}

if (verbose) {
		std::cout << "Masses:\n\tfrag1 = " << masses[0] << "\tfrag2 = " << masses[1] << "\n";
		std::cout << "\tRecoil = " << recoilMass << " + " << recoilEx << " = " << m_recoil << "\n";
		std::cout << "Decay constants:\n";
		std::cout << "\tEcm = " << expectedCMValues.Ecm << std::endl;
		std::cout << "\tVcm frag1 = " << expectedCMValues.vcm_frag1 << "\tKEcm frag1 = " << expectedCMValues.kecm_frag1 << std::endl;
		std::cout << "\tVcm frag2 = " << expectedCMValues.vcm_frag2 << "\tKEcm frag2 = " << expectedCMValues.kecm_frag2 << std::endl;
	}
}

int InvMass_Mult2::CountPermPasses(){
	int numPasses = 0;
	for(int i=0; i<2; i++) numPasses += caseResults[i].permPasses;
	return numPasses;
}

void InvMass_Mult2::ClearEventResults(){
	for(int i=0; i<2; i++) caseResults[i].Reset();
}
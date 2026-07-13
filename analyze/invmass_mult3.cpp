#include "invmass_mult3.h"
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

InvMass_Mult3::InvMass_Mult3()
	: outfile(nullptr), intermediateMass(0), recoilMass(0), recoilEx(0), beamEnergyMeV(0), intermediateEx(0){

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

	TString leaflist =  "SPSEnergy/D:SPSTheta/D:SPSPhi/D:SPS_Ex/D:" 
						"imIM/D:imEx/D:reconEx/D:imVCM/D:imKECM/D:imTHCM/D:imPHCM/D:imComp[3]/D:"
						"f1VCM/D:f1KECM/D:f1THCM/D:f1PHCM/D:f1Comp[3]/D:"
						"f2VCM/D:f2KECM/D:f2THCM/D:f2PHCM/D:f2Comp[3]/D:"
						"f3VCM/D:f3KECM/D:f3THCM/D:f3PHCM/D:f3Comp[3]/D:"
						"ecm1/D:ecm2/D:"
						"boost1[3]/D:boost2[3]/D:"
						"relLabAngle_intfrag1/D:relLabAngle_frag2frag3/D:"
						"IM2_int/D:"
						"exp_ecm1/D:exp_ecm2/D:"
						"exp_imVCM/D:exp_imKECM/D:"
						"exp_f1VCM/D:exp_f1KECM/D:"
						"exp_f2VCM/D:exp_f2KECM/D:"
						"exp_f3VCM/D:exp_f3KECM/D:"
						"permPasses/O";

	for(int i=0; i<6; i++){
		//create a branch for each permutation
		outtree->Branch("_"+permNames[i], &caseResults[i], leaflist);
	}

	for(auto &cn : permNames){
		TDirectory *permdir = outfile->mkdir(cn);
		//dir->cd();

		TDirectory *ungateDir = permdir->mkdir("ungated");
		TDirectory *gateDir = permdir->mkdir("gated");

		groups_ungated[cn] = new permHisto_mult3(cn + "_ungated", ungateDir);
		groups_gated[cn] = new permHisto_mult3(cn + "_gated", gateDir);

		outfile->cd();
	}

	hPermCounter = new TH1D("hPermCounter", "hPermCounter", 7, -0.5, 6.5);
	hPermCounter_gated = new TH1D("hPermCounter_gated", "hPermCounter_gated", 7, -0.5, 6.5);
	hSortedIntermediateExIMvsSPS = new TH2D("hSortedIntermediateExIMvsSPS", "Intermediate Ex (IM) vs Recoil Ex (SPS);SPS (MeV);IM (MeV)", 200, -1, 7, 200, -1, 7);
	hSortedDalitz = new TH2D("hSortedDalitz", "Dalitz m_{12}^{2} vs m_{23}^{2}", 1520/4, 55.55e6, 55.92e6, 1520, 0, 55.92e6);
	hSortedPermutations = new TH1D("hSortedPermutations", "Permutation Picks", 6, -0.5, 5.5);
	hSortedIMRecEx = new TH1D("hSortedIMRecEx", "Rec Ex (IM)", 300, -5, 7);
	hSortedIMRecEx_gate8Be = new TH1D("hSortedIMRecEx_gate8Be", "Rec Ex (IM) - Gate IM 8Be", 300, -5, 7);
	hSortedIMRecEx_gate5Li = new TH1D("hSortedIMRecEx_gate5Li", "Rec Ex (IM) - Gate IM 5Li", 300, -5, 7);
	hSABRESumE_vs_ExSPS = new TH2D("hSABRESumE_vs_ExSPS", "hSABRESumE_vs_ExSPS", 275, -1, 10, 500, 0, 10);
	const char* labels[6] = {"012","021","102","120","201","210"};
	for(int i=0; i<6; i++){
		hSortedPermutations->GetXaxis()->SetBinLabel(i+1, labels[i]);
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

	//calculate "above masses[0] threshold" in recoil here:
	
	recoilM0Threshold = intermediateMass + masses[0];// - recoilMass;

	//calculate "above masses[1] threshold" in intermediate here:
	intermediateM1Threshold = masses[1] + masses[2];// - intermediateMass;

	beamEnergyMeV = hypo.beamEnergyMeV;
	std::cout << "setting beamEnergyMeV to hypo.beamEnergyMeV (= " << hypo.beamEnergyMeV << ")...did it work?\n";

	//SetExpectedCMValues();
}

//Event assumes theta, phi in degrees and E in MeV
//update this to return reconstructed recoil excitation energy in MeV
std::array<double,6> InvMass_Mult3::AnalyzeEvent(double E[3], double theta[3], double phi[3], double SPSEnergy, double SPSTheta, double SPSPhi, double SPS_Ex, bool updateIntermediateEx){ 

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
			//double mom = std::sqrt(2*mass*E[hitindex]);//newtonian approximation
			double mom = std::sqrt(E[hitindex]*(E[hitindex]+2.*mass));//relativistic

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
		TLorentzVector frag1_boost2 = lv[0];//this is to be boosted into CM frame of second decay for use in calculating theta2h
		TLorentzVector frag2 = lv[1];
		TLorentzVector frag3 = lv[2];

		caseResults[permIndex].SPSEnergy = SPSEnergy;
		caseResults[permIndex].SPSTheta = SPSTheta;
		caseResults[permIndex].SPSPhi = SPSPhi;
		caseResults[permIndex].SPS_Ex = SPS_Ex;

		double PLabVect_beam_mag = std::sqrt(2.*hypothesis.mass_beam*beamEnergyMeV);
		TVector3 PLabVect_beam = {0., 0., PLabVect_beam_mag};

		double PLabVect_ejectile_mag = std::sqrt(2.*hypothesis.mass_ejectile*SPSEnergy);
		TVector3 PLabVect_ejectile = {PLabVect_ejectile_mag*std::sin(SPSTheta*DEGRAD)*std::cos(SPSPhi*DEGRAD),
									  PLabVect_ejectile_mag*std::sin(SPSTheta*DEGRAD)*std::sin(SPSPhi*DEGRAD),
									  PLabVect_ejectile_mag*std::cos(SPSTheta*DEGRAD)};

		TVector3 missingmomentum = PLabVect_beam - PLabVect_ejectile - (lv[0] + lv[1] + lv[2]).Vect();

		for(int i=0; i<3; i++) caseResults[permIndex].MissingMomentumComp[i] = missingmomentum[i];

		caseResults[permIndex].MissingMomentumMag = missingmomentum.Mag();

		caseResults[permIndex].PLabTotal_ResDecayParticles = (lv[0] + lv[1] + lv[2]).P();
		caseResults[permIndex].PLabTotal_Beam = std::sqrt(2*hypothesis.mass_beam*beamEnergyMeV);
		caseResults[permIndex].PLabTotal_Ejectile = std::sqrt(2*hypothesis.mass_ejectile*SPSEnergy);
		// std::cout << "PLabTotal_Beam = " << caseResults[permIndex].PLabTotal_Beam << "\n";
		// std::cout << "PLabTotal_ResDecayParticles = " << caseResults[permIndex].PLabTotal_ResDecayParticles << "\n";
		// std::cout << "PLabTotal_Ejectile = " << caseResults[permIndex].PLabTotal_Ejectile << "\n";
		// std::cout << "Missing momentum = " << missingmomentum.Mag() << "\n\n";
		//std::cout << "hypothesis.mass_beam = " << hypothesis.mass_beam << ", beamEnergyMeV = " << beamEnergyMeV << "\n";
		caseResults[permIndex].E2meas = E[indices[2]];//indices[2] is the hit index assigned to particle 2 (contents of indices[2] depends on current permutation)

		caseResults[permIndex].IM2_int = intermediate.M2();

		caseResults[permIndex].relLabAngle_intfrag1 = intermediate.Vect().Angle(frag1.Vect())*RADDEG;
		caseResults[permIndex].relLabAngle_frag2frag3 = frag2.Vect().Angle(frag3.Vect())*RADDEG;

		//calculate excitation energy:
		double intermediateEx = intermediate.M() - intermediateMass;
		double intermediateEnergyAboveM1Thresh = intermediate.M() - intermediateM1Threshold;
		double Ex = recoil.M() - recoilMass;
		double recoilEnergyAboveM0Thresh = recoil.M() - recoilM0Threshold;
		recoilExs[permIndex] = Ex;

		//let's check if we should update CM variables by using the IM-calcualted intermediate excitation energy:
		if(updateIntermediateEx){
			SetIntermediateEx(intermediateEx);
			//test for dynamic4:
			//SetRecoilEx(Ex);
		}

		caseResults[permIndex].permName = name;

		caseResults[permIndex].expected = expectedCMValues;

		caseResults[permIndex].intermediateIM = intermediate.M();
		caseResults[permIndex].intermediateEx = intermediateEx;
		caseResults[permIndex].reconEx = Ex;

		caseResults[permIndex].intermediateEnergyAboveM1Thresh = intermediateEnergyAboveM1Thresh;
		caseResults[permIndex].recoilEnergyAboveM0Thresh = recoilEnergyAboveM0Thresh;

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

		TVector3 intermediatevcm_vect = ((1/intermediate.Energy())*intermediate.Vect());
		double intermediatevcm = intermediatevcm_vect.Mag();
		//double intermediatevcm = intermediate.BoostVector().Mag();
		double intermediatekecm = 0.5*intermediateMass*intermediatevcm*intermediatevcm;
		double intermediatethetacm = RADDEG*std::acos(intermediate.Vect().Z()/intermediate.Vect().Mag());
		double intermediatephicm = RADDEG*std::atan2(intermediate.Vect().Y(), intermediate.Vect().X());
		if(intermediatephicm < 0) intermediatephicm += 360.;

		caseResults[permIndex].intermediatevcm = intermediatevcm;
		caseResults[permIndex].intermediatekecm = intermediatekecm;
		caseResults[permIndex].intermediatethetacm = intermediatethetacm;
		caseResults[permIndex].intermediatephicm = intermediatephicm;
		for(int i=0; i<3; i++) caseResults[permIndex].intermediateComp[i] = intermediatevcm_vect[i];

		TVector3 frag1vcm_vect = ((1/frag1.Energy())*frag1.Vect());
		double frag1vcm = frag1vcm_vect.Mag();
		//double frag1vcm = frag1.BoostVector().Mag();
		double frag1kecm = 0.5*masses[0]*frag1vcm*frag1vcm;
		double frag1thetacm = RADDEG*std::acos(frag1.Vect().Z()/frag1.Vect().Mag());
		double frag1phicm = RADDEG*std::atan2(frag1.Vect().Y(), frag1.Vect().X());
		if(frag1phicm < 0) frag1phicm += 360.;

		caseResults[permIndex].frag1vcm = frag1vcm;
		caseResults[permIndex].frag1kecm = frag1kecm;
		caseResults[permIndex].frag1thetacm = frag1thetacm;
		caseResults[permIndex].frag1phicm = frag1phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag1Comp[i] = frag1vcm_vect[i];

		//determine ecm1:
		double ecm1 = intermediatekecm + frag1kecm + intermediateEx;
		caseResults[permIndex].ecm1 = ecm1;


		//begin analysis of second decay step: intermediate -> frag2 + frag3
		//boost the lab-measured frag2 and frag3 into the frame of the intermediate:
		frag2.Boost(boost2);
		frag3.Boost(boost2);
		frag1_boost2.Boost(boost2);

		//calculate angle between boosted boosted frag2 and the boost vector of the intermediate:
		caseResults[permIndex].Theta2h = RADDEG*boost2.Angle(frag2.Vect());
		caseResults[permIndex].CosTheta2h = std::cos(boost2.Angle(frag2.Vect()));
		// double costheta2h_numerator = -2.*(frag2+frag3).M2()*(frag1_boost2+frag2).M2() + (frag2+frag3).M2()*(frag2.M2() + frag3.M2() - (frag2+frag3).M2()) + frag1_boost2.M2()*(frag2.M2() - frag3.M2() + (frag2+frag3).M2()) - (frag1_boost2+frag2+frag3).M2()*(frag2.M2() - frag3.M2() - (frag2+frag3).M2());
		// double costheta2h_denominator = std::sqrt((std::pow((frag1_boost2+frag2+frag3).M2(),2) + std::pow(frag1_boost2.M2() - (frag2+frag3).M2(),2) - 2.*(frag1_boost2+frag2+frag3).M2()*(frag1_boost2.M2()-(frag2+frag3).M2()))*(std::pow(frag2.M2(),2) + std::pow((frag2+frag3).M2() - frag3.M2(),2) - 2.*frag2.M2()*((frag2+frag3).M2() + frag3.M2())));
		// caseResults[permIndex].permTheta2h = costheta2h_numerator / costheta2h_denominator;

		TVector3 frag2vcm_vect = ((1/frag2.Energy())*frag2.Vect());
		double frag2vcm = frag2vcm_vect.Mag();
		//double frag2vcm = frag2.BoostVector().Mag();
		double frag2kecm = 0.5*masses[1]*frag2vcm*frag2vcm;
		double frag2thetacm = RADDEG*std::acos(frag2.Vect().Z()/frag2.Vect().Mag());
		double frag2phicm = RADDEG*std::atan2(frag2.Vect().Y(), frag2.Vect().X());
		if(frag2phicm < 0) frag2phicm += 360.;

		caseResults[permIndex].frag2vcm = frag2vcm;
		caseResults[permIndex].frag2kecm = frag2kecm;
		caseResults[permIndex].frag2thetacm = frag2thetacm;
		caseResults[permIndex].frag2phicm = frag2phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag2Comp[i] = frag2vcm_vect[i];

		TVector3 frag3vcm_vect = ((1/frag3.Energy())*frag3.Vect());
		double frag3vcm = frag3vcm_vect.Mag();
		//double frag3vcm = frag3.BoostVector().Mag();
		double frag3kecm = 0.5*masses[2]*frag3vcm*frag3vcm;
		double frag3thetacm = RADDEG*std::acos(frag3.Vect().Z()/frag3.Vect().Mag());
		double frag3phicm = RADDEG*std::atan2(frag3.Vect().Y(), frag3.Vect().X());
		if(frag3phicm < 0) frag3phicm += 360.;

		caseResults[permIndex].frag3vcm = frag3vcm;
		caseResults[permIndex].frag3kecm = frag3kecm;
		caseResults[permIndex].frag3thetacm = frag3thetacm;
		caseResults[permIndex].frag3phicm = frag3phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag3Comp[i] = frag3vcm_vect[i];

		//determine ecm2:
		double ecm2 = frag2kecm + frag3kecm;
		caseResults[permIndex].ecm2 = ecm2;

		caseResults[permIndex].exp_ecm1 = expectedCMValues.Ecm1;
		caseResults[permIndex].exp_ecm2 = expectedCMValues.Ecm2;
		caseResults[permIndex].exp_imVCM = expectedCMValues.vcm_intermediate;
		caseResults[permIndex].exp_imKECM = expectedCMValues.kecm_intermediate;
		caseResults[permIndex].exp_f1VCM  = expectedCMValues.vcm_frag1;
		caseResults[permIndex].exp_f1KECM = expectedCMValues.kecm_frag1;
		caseResults[permIndex].exp_f2VCM  = expectedCMValues.vcm_frag2;
		caseResults[permIndex].exp_f2KECM = expectedCMValues.kecm_frag2;
		caseResults[permIndex].exp_f3VCM  = expectedCMValues.vcm_frag3;
		caseResults[permIndex].exp_f3KECM = expectedCMValues.kecm_frag3;

		caseResults[permIndex].m12sq = (lv[0] + lv[1]).M2();
		caseResults[permIndex].m23sq = (lv[1] + lv[2]).M2();

		//increment permIndex here!
		permIndex += 1;

	}

	return recoilExs;
}

void InvMass_Mult3::FillTree(){
	if(outtree) outtree->Fill();
}

void InvMass_Mult3::FillEventHistograms(double SPS_Ex){
	for(int caseNum=0; caseNum<6; caseNum++){
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
		//fillAll("intermediateIM", res.intermediateIM);
		fillAll("intermediateEx", res.intermediateEx);
		fillAll("intermediateEnergyAboveM1Thresh", res.intermediateEnergyAboveM1Thresh);
		//std::cout << "intermediateEnergyAboveM1Thresh = " << res.intermediateEnergyAboveM1Thresh << "\n";
		fillAll("RecoilEx", res.reconEx);
		fillAll("RecoilEnergyAboveM0Thresh", res.recoilEnergyAboveM0Thresh);
		//std::cout << "recoilEnergyAboveM0Thresh = " << res.recoilEnergyAboveM0Thresh << "\n";
		fillAll2D("RecoilEx_IMvsSPS", SPS_Ex, res.reconEx);
		fillAll("RecoilExDif", SPS_Ex - res.reconEx);
		fillAll2D("intermediateExIMvsSPS", SPS_Ex, res.intermediateEx); 

		fillAll("MissingMomentum", res.MissingMomentumMag);
		//std::cout << "missing momentum = " << res.PLabTotal_Beam << " - " << res.PLabTotal_ResDecayParticles << " = " << res.PLabTotal_Beam - res.PLabTotal_ResDecayParticles << std::endl;

		//intermediate CM
		fillAll("intermediatevcm_meas", res.intermediatevcm);
		//fillAll("intermediatevcm_expect", res.expected.vcm_intermediate);
		//fillAll("intermediatevcm_delta", res.intermediatevcm - res.expected.vcm_intermediate);
		fillAll2D("intermediatevcm_TransverseVSLongitudinal", std::abs(res.intermediateComp[2]), std::sqrt(res.intermediateComp[0]*res.intermediateComp[0] + res.intermediateComp[1]*res.intermediateComp[1]));

		fillAll("intermediatekecm_meas", res.intermediatekecm);
		//fillAll("intermediatekecm_expect", res.expected.kecm_intermediate);
		//fillAll("intermediatekecm_delta", res.intermediatekecm - res.expected.kecm_intermediate);


		fillAll("intermediatethetacm", res.intermediatethetacm);
		fillAll("intermediatephicm", res.intermediatephicm);
		fillAll2D("intermediatethetacmvsphicm", res.intermediatephicm, res.intermediatethetacm);
		fillAll2D("intermediatevcmVSthetacm",res.intermediatethetacm, res.intermediatevcm);
		fillAll2D("intermediatekecmVSthetacm",res.intermediatethetacm, res.intermediatekecm);

		//fill frag1-3
		auto fillFrag = [&](int i, double vcm, double comps[3], double kecm, double theta, double phi, double expV, double expK){
			TString f = Form("frag%d",i);
			fillAll(f+"vcm_meas",vcm);
			//fillAll(f+"vcm_expect",expV);
			//fillAll(f+"vcm_delta",vcm-expV);
			fillAll2D(f+"vcm_TransverseVSLongitudinal", std::abs(comps[2]), std::sqrt(comps[0]*comps[0] + comps[1]*comps[1]));

			fillAll(f+"kecm_meas",kecm);
			//fillAll(f+"kecm_expect",expK);
			//fillAll(f+"kecm_delta",kecm-expK);

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
		//fillAll("ecm1_expect", res.expected.Ecm1);
		//fillAll("ecm1_delta", res.ecm1 - res.expected.Ecm1);
		fillAll2D("ecm1measVSintermediatethetacm", res.intermediatethetacm, res.ecm1);
		fillAll2D("ecm1measVSfrag1thetacm", res.frag1thetacm, res.ecm1);
		fillAll2D("ecm1measVSfrag2thetacm", res.frag2thetacm, res.ecm1);
		fillAll2D("ecm1measVSfrag3thetacm", res.frag3thetacm, res.ecm1);
		fillAll("decay1_VCM", std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1] + res.boost1[2]*res.boost1[2]));
		fillAll2D("decay1_VCM_TransverseVSLongitudinal", std::abs(res.boost1[2]), std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1]));
		fillAll("decay1_thetaCMsum", res.intermediatethetacm + res.frag1thetacm);
		fillAll("decay1_phiCMdiff", std::abs(res.intermediatephicm - res.frag1phicm));
		fillAll("decay1_relLabAngle", res.relLabAngle_intfrag1);
		fillAll("ecm2_meas", res.ecm2);
		//fillAll("ecm2_expect", res.expected.Ecm2);
		//fillAll("ecm2_delta", res.ecm2 - res.expected.Ecm2);
		fillAll2D("ecm2measVSintermediatethetacm", res.intermediatethetacm, res.ecm2);
		fillAll2D("ecm2measVSfrag1thetacm", res.frag1thetacm, res.ecm2);
		fillAll2D("ecm2measVSfrag2thetacm", res.frag2thetacm, res.ecm2);
		fillAll2D("ecm2measVSfrag3thetacm", res.frag3thetacm, res.ecm2);
		fillAll("decay2_VCM", std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1] + res.boost2[2]*res.boost2[2]));
		fillAll2D("decay2_VCM_TransverseVSLongitudinal", std::abs(res.boost2[2]), std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1]));
		fillAll("decay2_thetaCMsum", res.frag2thetacm + res.frag3thetacm);
		fillAll("decay2_phiCMdiff", std::abs(res.frag2phicm - res.frag3phicm));
		fillAll("decay2_relLabAngle", res.relLabAngle_frag2frag3);

		fillAll2D("decay2VSdecay1_relLabAngle", res.relLabAngle_intfrag1, res.relLabAngle_frag2frag3);

		fillAll2D("intermediatevcmVSfrag1vcm", res.frag1vcm, res.intermediatevcm);
		fillAll2D("intermediatekecmVSfrag1kecm", res.frag1kecm, res.intermediatekecm);
		fillAll2D("frag2vcmVSfrag3vcm", res.frag3vcm, res.frag2vcm);
		fillAll2D("frag2kecmVSfrag3kecm", res.frag3kecm, res.frag2kecm);
		fillAll2D("ecm1VSecm2", res.ecm2, res.ecm1);
		//fillAll2D("ecm1deltaVSecm2delta", res.ecm2 - res.expected.Ecm2, res.ecm1 - res.expected.Ecm1);

		//std::cout << "m12sq = " << res.m12sq << "\tm23sq = " << res.m23sq << std::endl;
		fillAll2D("dalitz_m12_vs_m23", res.m23sq, res.m12sq);

		fillAll("Theta2h", res.Theta2h);
		fillAll("CosTheta2h", res.CosTheta2h);
		fillAll2D("CosTheta2h_vs_m12sq", res.m12sq, res.CosTheta2h);
		fillAll2D("Theta2h_vs_DetE2", res.E2meas, res.Theta2h);
		fillAll2D("CosTheta2h_vs_DetE2", res.E2meas, res.CosTheta2h);
		fillAll2D("CosTheta2h_vs_MissingMomentum", res.MissingMomentumMag, res.CosTheta2h);
	}

	//if(outtree) outtree->Fill();
}

void InvMass_Mult3::FillGatedEventHistograms(double SPS_Ex){
	for(int caseNum=0; caseNum<6; caseNum++){
		TString cn = permNames.at(caseNum);
		auto& res = caseResults[caseNum];

		if(res.permName == "NONE") continue;

		// 1. Pack metrics and check if this specific permutation passes the gate
		std::map<TString, double> metrics1D = PackMetrics1D(res);
		std::map<TString, std::pair<double,double>> metrics2D = PackPoints2D(res);

		if (fCutHandler.PassAll1D(metrics1D) && fCutHandler.PassAll2D(metrics2D)) {

			// 2. Define lambdas to fill the gated histograms
			auto fillAllGated = [&](TString key, double x){
				groups_gated[cn]->Fill(key, x);
				groups_gated["allCases"]->Fill(key, x);
			};

			auto fillAll2DGated = [&](TString key, double x, double y){
				groups_gated[cn]->Fill(key, x, y);
				groups_gated["allCases"]->Fill(key, x, y);
			};

			// 3. Replicate the exact histogram filling logic from FillEventHistograms
			fillAllGated("intermediateEx", res.intermediateEx);
			fillAllGated("intermediateEnergyAboveM1Thresh", res.intermediateEnergyAboveM1Thresh);
			fillAllGated("RecoilEx", res.reconEx);
			fillAllGated("RecoilEnergyAboveM0Thresh", res.recoilEnergyAboveM0Thresh);
			fillAll2DGated("RecoilEx_IMvsSPS", SPS_Ex, res.reconEx);
			fillAllGated("RecoilExDif", SPS_Ex - res.reconEx);
			fillAll2DGated("intermediateExIMvsSPS", SPS_Ex, res.intermediateEx);

			fillAllGated("MissingMomentum", res.MissingMomentumMag);

			// Intermediate CM
			fillAllGated("intermediatevcm_meas", res.intermediatevcm);
			fillAll2DGated("intermediatevcm_TransverseVSLongitudinal", std::abs(res.intermediateComp[2]), std::sqrt(res.intermediateComp[0]*res.intermediateComp[0] + res.intermediateComp[1]*res.intermediateComp[1]));
			fillAllGated("intermediatekecm_meas", res.intermediatekecm);

			fillAllGated("intermediatethetacm", res.intermediatethetacm);
			fillAllGated("intermediatephicm", res.intermediatephicm);
			fillAll2DGated("intermediatethetacmvsphicm", res.intermediatephicm, res.intermediatethetacm);
			fillAll2DGated("intermediatevcmVSthetacm", res.intermediatethetacm, res.intermediatevcm);
			fillAll2DGated("intermediatekecmVSthetacm", res.intermediatethetacm, res.intermediatekecm);

			// Fragment helper lambda
			auto fillFragGated = [&](int i, double vcm, double comps[3], double kecm, double theta, double phi){
				TString f = Form("frag%d", i);
				fillAllGated(f+"vcm_meas", vcm);
				fillAll2DGated(f+"vcm_TransverseVSLongitudinal", std::abs(comps[2]), std::sqrt(comps[0]*comps[0] + comps[1]*comps[1]));
				fillAllGated(f+"kecm_meas", kecm);
				fillAllGated(f+"thetacm", theta);
				fillAllGated(f+"phicm", phi);
				fillAll2DGated(f+"thetacmvsphicm", phi, theta);
				fillAll2DGated(f+"vcmVSthetacm", theta, vcm);
				fillAll2DGated(f+"kecmVSthetacm", theta, kecm);
			};

			fillFragGated(1, res.frag1vcm, res.frag1Comp, res.frag1kecm, res.frag1thetacm, res.frag1phicm);
			fillFragGated(2, res.frag2vcm, res.frag2Comp, res.frag2kecm, res.frag2thetacm, res.frag2phicm);
			fillFragGated(3, res.frag3vcm, res.frag3Comp, res.frag3kecm, res.frag3thetacm, res.frag3phicm);

			// Decays
			fillAllGated("ecm1_meas", res.ecm1);
			fillAll2DGated("ecm1measVSintermediatethetacm", res.intermediatethetacm, res.ecm1);
			fillAll2DGated("ecm1measVSfrag1thetacm", res.frag1thetacm, res.ecm1);
			fillAll2DGated("ecm1measVSfrag2thetacm", res.frag2thetacm, res.ecm1);
			fillAll2DGated("ecm1measVSfrag3thetacm", res.frag3thetacm, res.ecm1);
			fillAllGated("decay1_VCM", std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1] + res.boost1[2]*res.boost1[2]));
			fillAll2DGated("decay1_VCM_TransverseVSLongitudinal", std::abs(res.boost1[2]), std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1]));
			fillAllGated("decay1_thetaCMsum", res.intermediatethetacm + res.frag1thetacm);
			fillAllGated("decay1_phiCMdiff", std::abs(res.intermediatephicm - res.frag1phicm));
			fillAllGated("decay1_relLabAngle", res.relLabAngle_intfrag1);

			fillAllGated("ecm2_meas", res.ecm2);
			fillAll2DGated("ecm2measVSintermediatethetacm", res.intermediatethetacm, res.ecm2);
			fillAll2DGated("ecm2measVSfrag1thetacm", res.frag1thetacm, res.ecm2);
			fillAll2DGated("ecm2measVSfrag2thetacm", res.frag2thetacm, res.ecm2);
			fillAll2DGated("ecm2measVSfrag3thetacm", res.frag3thetacm, res.ecm2);
			fillAllGated("decay2_VCM", std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1] + res.boost2[2]*res.boost2[2]));
			fillAll2DGated("decay2_VCM_TransverseVSLongitudinal", std::abs(res.boost2[2]), std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1]));
			fillAllGated("decay2_thetaCMsum", res.frag2thetacm + res.frag3thetacm);
			fillAllGated("decay2_phiCMdiff", std::abs(res.frag2phicm - res.frag3phicm));
			fillAllGated("decay2_relLabAngle", res.relLabAngle_frag2frag3);

			fillAll2DGated("decay2VSdecay1_relLabAngle", res.relLabAngle_intfrag1, res.relLabAngle_frag2frag3);
			fillAll2DGated("intermediatevcmVSfrag1vcm", res.frag1vcm, res.intermediatevcm);
			fillAll2DGated("intermediatekecmVSfrag1kecm", res.frag1kecm, res.intermediatekecm);
			fillAll2DGated("frag2vcmVSfrag3vcm", res.frag3vcm, res.frag2vcm);
			fillAll2DGated("frag2kecmVSfrag3kecm", res.frag3kecm, res.frag2kecm);
			fillAll2DGated("ecm1VSecm2", res.ecm2, res.ecm1);

			fillAll2DGated("dalitz_m12_vs_m23", res.m23sq, res.m12sq);

			fillAllGated("Theta2h", res.Theta2h);
			fillAllGated("CosTheta2h", res.CosTheta2h);
			fillAll2DGated("CosTheta2h_vs_m12sq", res.m12sq, res.CosTheta2h);
			fillAll2DGated("Theta2h_vs_DetE2", res.E2meas, res.Theta2h);
			fillAll2DGated("CosTheta2h_vs_DetE2", res.E2meas, res.CosTheta2h);
			fillAll2DGated("CosTheta2h_vs_MissingMomentum", res.MissingMomentumMag, res.CosTheta2h);
		}
	}
}

void InvMass_Mult3::FillPermCounter(bool gated){
	if(gated){
		hPermCounter_gated->Fill(CountPermPasses());
	} else {
		hPermCounter->Fill(CountPermPasses());
	}
}

void InvMass_Mult3::FillSortedHisto(double SPS_Ex){
	//set new dummy array equal to caseResults and sort it by difference between IM Ex and SPS Ex
	Results sorted_caseResults[6];
	for(int i=0; i<6; i++) sorted_caseResults[i] = caseResults[i];

	std::sort(std::begin(sorted_caseResults), std::end(sorted_caseResults), [SPS_Ex](const Results& a, const Results& b) { 
		double diffA = std::abs(SPS_Ex - a.reconEx);
		double diffB = std::abs(SPS_Ex - b.reconEx);
		return diffA < diffB;
	});

	hSortedIntermediateExIMvsSPS->Fill(SPS_Ex, sorted_caseResults[0].intermediateEx);
	hSortedDalitz->Fill(sorted_caseResults[0].m23sq, sorted_caseResults[0].m12sq);
	hSortedIMRecEx->Fill(sorted_caseResults[0].reconEx);
	if(sorted_caseResults[0].intermediateEx >= -0.04 && sorted_caseResults[0].intermediateEx <= 0.04){
		hSortedIMRecEx_gate8Be->Fill(sorted_caseResults[0].reconEx);
	}
	if(sorted_caseResults[0].intermediateEx >= -0.8 && sorted_caseResults[0].intermediateEx <= 0.8){
		hSortedIMRecEx_gate5Li->Fill(sorted_caseResults[0].reconEx);
	}

	int permIndex = -1;
	TString perm = sorted_caseResults[0].permName;
	if(perm == "012") permIndex = 0;
	else if(perm == "021") permIndex = 1;
	else if(perm == "102") permIndex = 2;
	else if(perm == "120") permIndex = 3;
	else if(perm == "201") permIndex = 4;
	else if(perm == "210") permIndex = 5;

	hSortedPermutations->Fill(permIndex);
}

void InvMass_Mult3::FillSABRESumEVsSPS_Ex(double SPS_Ex, double SABREsumE){
	hSABRESumE_vs_ExSPS->Fill(SPS_Ex, SABREsumE);
}

void InvMass_Mult3::CloseAndWrite(){
	outfile->cd();
	// hPermCounter->Write();
	// hPermCounter_gated->Write();
	// hSortedIntermediateExIMvsSPS->Write();
	// hSortedDalitz->Write();
	// hSortedPermutations->Write();
	// hSortedIMRecEx->Write();
	// hSortedIMRecEx_gate8Be->Write();
	// hSortedIMRecEx_gate5Li->Write();
	// hSABRESumE_vs_ExSPS->Write();
	if(outfile && outfile->IsOpen()){
		outfile->Write();
		outfile->Close();
	}
	// delete hPermCounter;
	// delete hPermCounter_gated;
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

std::map<TString, double> InvMass_Mult3::PackMetrics1D(const Results& res) const {
	std::map<TString, double> metrics;

	metrics["intermediateEx"] = res.intermediateEx;
	metrics["reconEx"] = res.reconEx;
	metrics["intermediateIM"] = res.intermediateIM;
	metrics["intermediateEnergyAboveM1Thresh"] = res.intermediateEnergyAboveM1Thresh;
	metrics["recoilEnergyAboveM0Thresh"] = res.recoilEnergyAboveM0Thresh;
	metrics["E2meas"] = res.E2meas;

	metrics["intermediatevcm"] = res.intermediatevcm;
	metrics["intermediatekecm"] = res.intermediatekecm;
	metrics["frag1vcm"] = res.frag1vcm;
	metrics["frag1kecm"] = res.frag1kecm;
	metrics["frag2vcm"] = res.frag2vcm;
	metrics["frag2kecm"] = res.frag2kecm;
	metrics["frag3vcm"] = res.frag3vcm;
	metrics["frag3kecm"] = res.frag3kecm;

	metrics["Theta2h"] = res.Theta2h;
	metrics["CosTheta2h"] = res.CosTheta2h;
	metrics["relLabAngle_intfrag1"] = res.relLabAngle_intfrag1;
	metrics["relLabAngle_frag2frag3"] = res.relLabAngle_frag2frag3;

	return metrics;
}

std::map<TString, std::pair<double, double>> InvMass_Mult3::PackPoints2D(const Results& res) const {
	std::map<TString, std::pair<double,double>> points;

	points["RecoilEx_IMvsSPS"] = { res.SPS_Ex, res.reconEx };
	points["intermediateExIMvsSPS"] = { res.SPS_Ex, res.intermediateEx };
	points["intermediatevcm_TransverseVSLongitudinal"] = { std::abs(res.intermediateComp[2]), std::sqrt(res.intermediateComp[0]*res.intermediateComp[0] + res.intermediateComp[1]*res.intermediateComp[1])};
	points["frag1vcm_TransverseVSLongitudinal"] = { std::abs(res.frag1Comp[2]), std::sqrt(res.frag1Comp[0]*res.frag1Comp[0] + res.frag1Comp[1]*res.frag1Comp[1]) };
	points["frag2vcm_TransverseVSLongitudinal"] = { std::abs(res.frag2Comp[2]), std::sqrt(res.frag2Comp[0]*res.frag2Comp[0] + res.frag2Comp[1]*res.frag2Comp[1]) };
	points["frag3vcm_TransverseVSLongitudinal"] = { std::abs(res.frag3Comp[2]), std::sqrt(res.frag3Comp[0]*res.frag3Comp[0] + res.frag3Comp[1]*res.frag3Comp[1]) };
	points["CosTheta2h_vs_MissingMomentum"] = { res.MissingMomentumMag, res.CosTheta2h };
	points["ecm1VSecm2"] = { res.ecm2, res.ecm1 };
	points["dalitz_m12_vs_m23"] = { res.m23sq, res.m12sq };

	return points;
}

int InvMass_Mult3::CountPermPasses(){
	int passCount = 0;

	for(int i=0; i<6; i++){
		auto& res = caseResults[i];
		if(res.permName == "NONE"){
			res.permPasses = false;
			continue;
		}

		std::map<TString, double> metrics1D = PackMetrics1D(res);
		std::map<TString, std::pair<double,double>> metrics2D = PackPoints2D(res);

		if(fCutHandler.PassAll1D(metrics1D) && fCutHandler.PassAll2D(metrics2D)){
			res.permPasses = true;
			passCount++;
		} else {
			res.permPasses = false;
		}
	}

	return passCount;
}

void InvMass_Mult3::ClearEventResults(){
	for(int i=0; i<6; i++){
		caseResults[i].Reset();
	}
}
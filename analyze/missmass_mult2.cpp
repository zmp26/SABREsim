#include "missmass_mult2.h"
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>

	/*

		README BEFORE READING CODE!

		The purpose of this code is to perform kinematically constrained
		missing mass analysis of events in which only N-1 resonance decay
		particles are detected in coincidence with the reaction ejectile.
		This code specifically handles events from decays resulting in
		N=3 final particles of which only N-1=2 are detected in addition
		to the reaction ejectile.

		Assume some nuclear reaction A(a,b)B* in which an unbound resonance
		B* is populated and then decays via B*->C*+c where C* is another
		unbound resonance (may or may not be an excited state) and c is
		a stable resonance-decay particle. Shortly after being produced, C*
		undergoes a decay C*->d+e. Thus, the final state of the reaction 
		consists only of the outgoing reaction ejectile b and the stable
		resonance decay particles c,d,e.

		In the code, c,d,e are referred to as frag1, frag2, and frag3.
		So, we have:

		B* -> frag1 + intermediate
								  \
					          	   frag2 + frag3

		This code ALWAYS assumes that frag2+frag3 = intermediate, and
		further that intermeditae+frag1 = frag1+frag2+frag3 = recoil.

		Since this code analyzes events where we have only two frags
		and the ejectile, we have three main cases to check:

			I:		Measure frag2 + frag3	-->	Reconstruct frag1
			II:		Measure frag1 + frag2	--> Reconstruct frag3
			III:	Measure frag1 + frag3	--> Reconstruct frag2

		The code assigns the two EThetaPhi values passed in AnalyzeEvent
		to the fragments shown above. Because this code assumes no particle
		ID for resonance-decay particles, there two sub cases for the two
		possible assignments. these subcases are shown below:

			I:		Measure frag2 + frag3	-->	Reconstruct frag1
					Ia:		EThetaPhi[0] = frag2, EThetaPhi[1] = frag3
					Ib: 	EThetaPhi[0] = frag3, EThetaPhi[0] = frag2

			II:		Measure frag1 + frag2	--> Reconstruct frag3
					IIa:	EThetaPhi[0] = frag1, EThetaPhi[1] = frag2
					IIb:	EThetaPhi[0] = frag2, EThetaPhi[1] = frag1

			III:	Measure frag1 + frag3	--> Reconstruct frag2
					IIIa:	EThetaPhi[0] = frag1, EThetaPhi[1] = frag3
					IIIb:	EThetaPhi[0] = frag3, EThetaPhi[1] = frag1

		The "permutations" in permNames and pMap refer to these 6 individual
		cases, which are "permutations" of the possible ways to assign our two
		EThetaPhi values to any subset of two paticles from the decay hypothesis.
		These labels are telling you which particle is assumed to be missing (and
		thus reconstructed), and the a/b tells you which of the two possible sub
		cases. See table here:

			recon_f1a		-->			case Ia
			recon_f1b		-->			case Ib

			recon_f2a		-->			case IIa
			recon_f2b		-->			case IIb

			recon_f3a		-->			case IIIa
			recon_f3b		-->			case IIIb

			allCases		-->			simply just the above 6 summed on one histogram

		This code tries all 6 cases for every event passed to AnalyzeEvent, and
		then plots histograms of the results in the assigned TDirectory.

	*/

MissMass_Mult2::MissMass_Mult2()
	: outfile(nullptr), beamMass(0), targetMass(0), ejectileMass(0), beamEnergyMeV(0), intermediateMass(0), recoilMass(0), intermediateEx(0){

	permNames = {"recon_f1a", "recon_f1b", "recon_f2a", "recon_f2b", "recon_f3a", "recon_f3b", "allCases"};

	pMap = {
		{permNames[0],{-1,0,1}},	//Ia
		{permNames[1],{-1,1,0}},	//Ib
		{permNames[2],{0,1,-1}},	//IIa
		{permNames[3],{1,0,-1}},	//IIb
		{permNames[4],{0,-1,1}},	//IIIa
		{permNames[5],{1,-1,0}}		//IIIb
	};

}

MissMass_Mult2::~MissMass_Mult2(){
	if(outfile && !outfile->IsZombie()){
		outfile->Close();
		delete outfile;
	}
}

void MissMass_Mult2::Init(const char* output_filename){

	outfile = new TFile(output_filename, "RECREATE");
	outtree = new TTree("MissMass_Mult2", "MissMass_Mult2");

	TString leaflist = "imIM/D:imEx/D:reconEx/D:imVCM/D:imKECM/D:imTHCM/D:imPHCM/D:imComp[3]/D:"
					   "f1VCM/D:f1KECM/D:f1THCM/D:f1PHCM/D:f1Comp[3]/D:"
					   "f2VCM/D:f2KECM/D:f2THCM/D:f2PHCM/D:f2Comp[3]/D:"
					   "f3VCM/D:f3KECM/D:f3THCM/D:f3PHCM/D:f3Comp[3]/D:"
					   "missingmass/D:"
					   "ecm1/D:ecm2/D:"
					   "boost1[3]/D:boost2[3]/D:"
					   "relLabAngle_intfrag1/D:relLabAngle_frag2frag3/D:"
					   //"IM2_01/D:IM2_12/D:IM2_20/D:"
					   "IM2_int/D:"
					   "exp_ecm1/D:exp_ecm2/D:"
					   "exp_imVCM/D:exp_imKECM/D:"
					   "exp_f1VCM/D:exp_f1KECM/D:"
					   "exp_f2VCM/D:exp_f2KECM/D:"
					   "exp_f3VCM/D:exp_f3KECM/D:"
					   "permPasses/O";

	for(int i=0; i<6; i++){
		outtree->Branch(permNames[i], &caseResults[i], leaflist);
	}

	for(auto& pn : permNames){
		TDirectory *permdir = outfile->mkdir(pn);

		TDirectory *ungateDir = permdir->mkdir("ungated");
		TDirectory *gateDir = permdir->mkdir("gated");

		groups_ungated[pn] = new permHistoMM_mult2(pn+"_ungated", ungateDir);
		groups_gated[pn] = new permHistoMM_mult2(pn+"_gated", gateDir);

		outfile->cd();
	}

	hPermCounter = new TH1D("hPermCounter", "hPermCounter", 7, -0.5, 6.5);
	hPermCounter_gated = new TH1D("hPermCounter_gated", "hPermCounter_gated", 7, -0.5, 6.5);
	hSortedIntermediateExIMvsSPS = new TH2D("hSortedIntermediateExIMvsSPS", "Intermediate Ex (IM) vs Recoil Ex (SPS);SPS (MeV);IM (MeV)", 200, -1, 7, 200, -1, 7);
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

void MissMass_Mult2::SetHypothesis(const Hypothesis4MM& hypo){
	hypothesis = hypo;

	for(int i=0; i<3; i++) masses[i] = hypo.masses[i];

	recoilMass = hypo.mass_recoil;
	intermediateMass = hypo.mass_intermediate;
	beamMass = hypo.mass_beam;
	targetMass = hypo.mass_target;
	ejectileMass = hypo.mass_ejectile;
	beamEnergyMeV = hypo.beamEnergyMeV;

}

std::array<double,6> MissMass_Mult2::AnalyzeEvent(double E[2], double theta[2], double phi[2], double SPSE, double SPSTheta, double SPSPhi, bool updateIntermediateEx){

	ClearEventResults();

	std::array<double,6> recoilExs;

	int permIndex = 0;
	for(auto const& [name, p] : pMap){

		TLorentzVector lv[3];//lv[0] = frag1, lv[1] = frag2, lv[2] = frag3
		int indices[3] = {p.i, p.j, p.k};

		//reconstruct missing particle from measured particles and ejectile
		TLorentzVector beam, target, ejectile, recoil;

		double beam_mom = std::sqrt(2*beamMass*beamEnergyMeV);
		beam.SetPxPyPzE(
				0,
				0,
				beam_mom,
				beamEnergyMeV+beamMass
			);

		target.SetPxPyPzE(
				0,
				0,
				0,
				targetMass
			);

		double ej_mom = std::sqrt(2*ejectileMass*SPSE);
		ejectile.SetPxPyPzE(
				ej_mom*std::sin(SPSTheta*DEGRAD)*std::cos(SPSPhi*DEGRAD),
				ej_mom*std::sin(SPSTheta*DEGRAD)*std::sin(SPSPhi*DEGRAD),
				ej_mom*std::cos(SPSTheta*DEGRAD),
				SPSE+ejectileMass
			);

		//recoil is simply initial - ejectile:
		recoil = beam + target - ejectile;

		//depending on which permIndex we are on, we assign the EThetaPhi to different particle assumptions (see note at top of file)
		double missingmass = -666.;
		switch(permIndex){
			double mass, mom;
			case 0: 
				//Ia, so assign [0] to frag2, [1] to frag3
				mass = masses[1];
				mom = std::sqrt(2*mass*E[0]);
				lv[1].SetPxPyPzE(
						mom*std::sin(theta[0]*DEGRAD)*std::cos(phi[0]*DEGRAD),
						mom*std::sin(theta[0]*DEGRAD)*std::sin(phi[0]*DEGRAD),
						mom*std::cos(theta[0]*DEGRAD),
						E[0] + mass
					);

				mass = masses[2];
				mom = std::sqrt(2*mass*E[1]);
				lv[2].SetPxPyPzE(
						mom*std::sin(theta[1]*DEGRAD)*std::cos(phi[1]*DEGRAD),
						mom*std::sin(theta[1]*DEGRAD)*std::sin(phi[1]*DEGRAD),
						mom*std::cos(theta[1]*DEGRAD),
						E[1] + mass
					);
				break;
			case 1:
				//Ib, so assign [0] to frag3, [1] to frag2
				mass = masses[2];
				mom = std::sqrt(2*mass*E[0]);
				lv[1].SetPxPyPzE(
						mom*std::sin(theta[0]*DEGRAD)*std::cos(phi[0]*DEGRAD),
						mom*std::sin(theta[0]*DEGRAD)*std::sin(phi[0]*DEGRAD),
						mom*std::cos(theta[0]*DEGRAD),
						E[0] + mass
					);

				mass = masses[1];
				mom = std::sqrt(2*mass*E[1]);
				lv[2].SetPxPyPzE(
						mom*std::sin(theta[1]*DEGRAD)*std::cos(phi[1]*DEGRAD),
						mom*std::sin(theta[1]*DEGRAD)*std::sin(phi[1]*DEGRAD),
						mom*std::cos(theta[1]*DEGRAD),
						E[1] + mass
					);
				break;
			case 2:
				//IIa, so assign [0] to frag1, [1] to frag2
				mass = masses[0];
				mom = std::sqrt(2*mass*E[0]);
				lv[0].SetPxPyPzE(
						mom*std::sin(theta[0]*DEGRAD)*std::cos(phi[0]*DEGRAD),
						mom*std::sin(theta[0]*DEGRAD)*std::sin(phi[0]*DEGRAD),
						mom*std::cos(theta[0]*DEGRAD),
						E[0] + mass
					);

				mass = masses[1];
				mom = std::sqrt(2*mass*E[1]);
				lv[1].SetPxPyPzE(
						mom*std::sin(theta[1]*DEGRAD)*std::cos(phi[1]*DEGRAD),
						mom*std::sin(theta[1]*DEGRAD)*std::sin(phi[1]*DEGRAD),
						mom*std::cos(theta[1]*DEGRAD),
						E[1] + mass
					);
				break;
			case 3:
				//IIb, so assign [0] to frag2, [1] to frag1
				mass = masses[1];
				mom = std::sqrt(2*mass*E[0]);
				lv[0].SetPxPyPzE(
						mom*std::sin(theta[0]*DEGRAD)*std::cos(phi[0]*DEGRAD),
						mom*std::sin(theta[0]*DEGRAD)*std::sin(phi[0]*DEGRAD),
						mom*std::cos(theta[0]*DEGRAD),
						E[0] + mass
					);

				mass = masses[0];
				mom = std::sqrt(2*mass*E[1]);
				lv[1].SetPxPyPzE(
						mom*std::sin(theta[1]*DEGRAD)*std::cos(phi[1]*DEGRAD),
						mom*std::sin(theta[1]*DEGRAD)*std::sin(phi[1]*DEGRAD),
						mom*std::cos(theta[1]*DEGRAD),
						E[1] + mass
					);
				break;
			case 4:
				//IIIa, so assign [0] to frag1, [1] to frag3
				mass = masses[0];
				mom = std::sqrt(2*mass*E[0]);
				lv[0].SetPxPyPzE(
						mom*std::sin(theta[0]*DEGRAD)*std::cos(phi[0]*DEGRAD),
						mom*std::sin(theta[0]*DEGRAD)*std::sin(phi[0]*DEGRAD),
						mom*std::cos(theta[0]*DEGRAD),
						E[0] + mass
					);

				mass = masses[2];
				mom = std::sqrt(2*mass*E[1]);
				lv[2].SetPxPyPzE(
						mom*std::sin(theta[1]*DEGRAD)*std::cos(phi[1]*DEGRAD),
						mom*std::sin(theta[1]*DEGRAD)*std::sin(phi[1]*DEGRAD),
						mom*std::cos(theta[1]*DEGRAD),
						E[1] + mass
					);
				break;
			case 5:
				//IIIb, so assign [0] to frag3, [1] to frag1
				mass = masses[2];
				mom = std::sqrt(2*mass*E[0]);
				lv[0].SetPxPyPzE(
						mom*std::sin(theta[0]*DEGRAD)*std::cos(phi[0]*DEGRAD),
						mom*std::sin(theta[0]*DEGRAD)*std::sin(phi[0]*DEGRAD),
						mom*std::cos(theta[0]*DEGRAD),
						E[0] + mass
					);

				mass = masses[0];
				mom = std::sqrt(2*mass*E[1]);
				lv[2].SetPxPyPzE(
						mom*std::sin(theta[1]*DEGRAD)*std::cos(phi[1]*DEGRAD),
						mom*std::sin(theta[1]*DEGRAD)*std::sin(phi[1]*DEGRAD),
						mom*std::cos(theta[1]*DEGRAD),
						E[1] + mass
					);
				break;
			default: break;
		}

		//no matter which permIndex we are in, the missing particle is always calculated the same way
		//first, we generally identify which one is missing:
		TLorentzVector lv_measuredSum;
		int missing_idx = -1;

		for(int m=0; m<3; m++){
			if(lv[m].E() > 0) lv_measuredSum += lv[m];
			else missing_idx = m;
		}

		caseResults[permIndex].SABREsumE = lv_measuredSum.E();

		if(missing_idx != -1) lv[missing_idx] = recoil - lv_measuredSum;
		missingmass = lv[missing_idx].M();

		TLorentzVector frag1, frag2, frag3;//just for reading simplicity instead of lv[#]
		frag1 = lv[0];
		frag2 = lv[1];
		frag3 = lv[2];

		TLorentzVector intermediate = frag2 + frag3;
		double intermediateIM = intermediate.M();

		double imEx = intermediateIM - intermediateMass;
		TLorentzVector lv_recoil_reconstructed = frag1 + intermediate;
		double reconEx = lv_recoil_reconstructed.M() - recoilMass;

		caseResults[permIndex].intermediateIM = intermediateIM;
		caseResults[permIndex].IM2_int = intermediate.M2();
		caseResults[permIndex].intermediateEx = imEx;
		caseResults[permIndex].reconEx = reconEx;

		caseResults[permIndex].relLabAngle_intfrag1 = intermediate.Vect().Angle(frag1.Vect());
		caseResults[permIndex].relLabAngle_frag2frag3 = frag2.Vect().Angle(frag3.Vect());

		if(updateIntermediateEx){
			SetIntermediateEx(imEx);
		}

		caseResults[permIndex].permName = name;
		caseResults[permIndex].expected = expectedCMValues;

		TVector3 boost1 = -recoil.BoostVector();
		TVector3 boost2 = -intermediate.BoostVector();

		caseResults[permIndex].boost1[0] = boost1.X();
		caseResults[permIndex].boost1[1] = boost1.Y();
		caseResults[permIndex].boost1[2] = boost1.Z();

		caseResults[permIndex].boost2[0] = boost2.X();
		caseResults[permIndex].boost2[1] = boost2.Y();
		caseResults[permIndex].boost2[2] = boost2.Z();

		//begin analysis of first decay step here: recoil->frag1+intermediate
		//boost the lab-measured intermediate and frag1 into the frame of the recoil:
		intermediate.Boost(boost1);
		frag1.Boost(boost1);

		double intermediatevcm = ((1/intermediate.Energy())*intermediate.Vect()).Mag();
		//double intermediatevcm = intermediate.BoostVector().Mag();
		double intermediatekecm = 0.5*intermediateMass*intermediatevcm*intermediatevcm;
		// double intermediatethetacm = RADDEG*std::acos(intermediate.Vect().Z()/intermediate.Vect().Mag());
		// double intermediatephicm = RADDEG*std::atan2(intermediate.Vect().Y(), intermediate.Vect().X());
		double intermediatethetacm = intermediate.Theta()*RADDEG;
		double intermediatephicm = intermediate.Phi()*RADDEG;
		if(intermediatephicm < 0) intermediatephicm += 360.;

		caseResults[permIndex].intermediatevcm = intermediatevcm;
		caseResults[permIndex].intermediatekecm = intermediatekecm;
		caseResults[permIndex].intermediatethetacm = intermediatethetacm;
		caseResults[permIndex].intermediatephicm = intermediatephicm;
		for(int i=0; i<3; i++) caseResults[permIndex].intermediateComp[i] = (-intermediate.BoostVector())(i);

		double frag1vcm = ((1/frag1.Energy())*frag1.Vect()).Mag();
		//double frag1vcm = frag1.BoostVector().Mag();
		double frag1kecm = 0.5*masses[0]*frag1vcm*frag1vcm;
		// double frag1thetacm = RADDEG*std::acos(frag1.Vect().Z()/frag1.Vect().Mag());
		// double frag1phicm = RADDEG*std::atan2(frag1.Vect().Y(), frag1.Vect().X());
		double frag1thetacm = frag1.Theta()*RADDEG;
		double frag1phicm = frag1.Phi()*RADDEG;
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
		//double frag2thetacm = RADDEG*std::acos(frag2.Vect().Z()/frag2.Vect().Mag());
		//double frag2phicm = RADDEG*std::atan2(frag2.Vect().Y(), frag2.Vect().X());
		double frag2thetacm = frag2.Theta()*RADDEG;
		double frag2phicm = frag2.Phi()*RADDEG;
		if(frag2phicm < 0) frag2phicm += 360.;

		caseResults[permIndex].frag2vcm = frag2vcm;
		caseResults[permIndex].frag2kecm = frag2kecm;
		caseResults[permIndex].frag2thetacm = frag2thetacm;
		caseResults[permIndex].frag2phicm = frag2phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag2Comp[i] = (-frag2.BoostVector())(i);

		double frag3vcm = ((1/frag3.Energy())*frag3.Vect()).Mag();
		//double frag3vcm = frag3.BoostVector().Mag();
		double frag3kecm = 0.5*masses[2]*frag3vcm*frag3vcm;
		//double frag3thetacm = RADDEG*std::acos(frag3.Vect().Z()/frag3.Vect().Mag());
		//double frag3phicm = RADDEG*std::atan2(frag3.Vect().Y(), frag3.Vect().X());
		double frag3thetacm = frag3.Theta()*RADDEG;
		double frag3phicm = frag3.Phi()*RADDEG;
		if(frag3phicm < 0) frag3phicm += 360.;

		caseResults[permIndex].frag3vcm = frag3vcm;
		caseResults[permIndex].frag3kecm = frag3kecm;
		caseResults[permIndex].frag3thetacm = frag3thetacm;
		caseResults[permIndex].frag3phicm = frag3phicm;
		for(int i=0; i<3; i++) caseResults[permIndex].frag3Comp[i] = (-frag3.BoostVector())(i);

		caseResults[permIndex].missingmass = missingmass;

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

		//increment permIndex:
		permIndex += 1;

		//temp debugging:
		// std::cout << "f1 ThetaCM=" << caseResults[permIndex].frag1thetacm <<", PhiCM = " << caseResults[permIndex].frag1phicm << std::endl;
		// std::cout << "f2 ThetaCM=" << caseResults[permIndex].frag2thetacm <<", PhiCM = " << caseResults[permIndex].frag2phicm << std::endl;
		// std::cout << "f3 ThetaCM=" << caseResults[permIndex].frag3thetacm <<", PhiCM = " << caseResults[permIndex].frag3phicm << std::endl;
		// std::cout << "f1 ThetaCM=" << frag1thetacm <<", PhiCM = " << frag1phicm << std::endl;
		// std::cout << "f2 ThetaCM=" << frag2thetacm <<", PhiCM = " << frag2phicm << std::endl;
		// std::cout << "f3 ThetaCM=" << frag3thetacm <<", PhiCM = " << frag3phicm << std::endl;
		// std::cout << "\n";

	}

	for(int i=0; i<6; i++){

		bool gate1 = false, gate2 = false, gate3 = false;

		auto& res = caseResults[i];

		switch(gate1_index){
			case NOCHECK:
				gate1 = true;
				break;
			case INTEXCHECK:
				//gate1 = intermediateExCheck;
				//gate1 = (res.intermediateEx >= gate1minmax.first && res.intermediateEx <= gate1minmax.second);
				gate1 = CheckGate1(res.intermediateEx);
				break;
			case INTVCMCHECK:
				//gate1 = intermediateVcmCheck;
				//gate1 = (res.intermediatevcm >= gate1minmax.first && res.intermediatevcm <= gate1minmax.second);
				gate1 = CheckGate1(res.intermediatevcm);
				break;
			case FRAG1VCMCHECK:
				//gate1 = (res.frag1vcm >= gate1minmax.first && res.frag1vcm <= gate1minmax.second);
				gate1 = CheckGate1(res.frag1vcm);
				break;
			default:
				gate1 = false;
				break;
		}

		switch(gate2_index){
			case NOCHECK:
				gate2 = true;
				break;
			case INTEXCHECK:
				//gate2 = intermediateExCheck;
				//gate2 = (res.intermediateEx >= gate2minmax.first && res.intermediateEx <= gate2minmax.second);
				gate2 = CheckGate2(res.intermediateEx);
				break;
			case INTVCMCHECK:
				//gate2 = intermediateVcmCheck;
				//gate2 = (res.intermediatevcm >= gate2minmax.first && res.intermediatevcm <= gate2minmax.second);
				gate2 = CheckGate2(res.intermediatevcm);
				break;
			case FRAG1VCMCHECK:
				//gate2 = frag1VcmCheck;
				//gate2 = (res.frag1vcm >= gate2minmax.first && res.frag1vcm <= gate2minmax.second);
				gate2 = CheckGate2(res.frag1vcm);
				break;
			default:
				gate2 = false;
				break;
		}

		if(gate1 && gate2) res.permPasses = true;
	}

	return recoilExs;

}

void MissMass_Mult2::FillTree(){
	if(outtree) outtree->Fill();
}

void MissMass_Mult2::FillEventHistograms(double SPS_Ex){
	for(int i=0; i<6; i++) FillSelectCaseHistograms(i, SPS_Ex);
}

void MissMass_Mult2::FillSelectCaseHistograms(int caseNum, double SPS_Ex){
	if(caseNum < 0 || caseNum > 5){
		std::cerr << "caseNum out of range (" << caseNum << ")\n";
		return;
	}

	TString pn = permNames.at(caseNum);
	auto& res = caseResults[caseNum];

	auto fillAll = [&](TString key, double x){
		groups_ungated[pn]->Fill(key,x);
		groups_ungated["allCases"]->Fill(key,x);
	};

	auto fillAll2D = [&](TString key, double x, double y){
		groups_ungated[pn]->Fill(key, x, y);
		groups_ungated["allCases"]->Fill(key, x, y);
	};

	//invariant mass and excitation energy histograms:
	fillAll("intermediateIM", res.intermediateIM);
	fillAll("intermediateEx", res.intermediateEx);
	fillAll("RecoilEx", res.reconEx);
	fillAll2D("RecoilEx_IMvsSPS", SPS_Ex, res.reconEx);
	fillAll("RecoilExDif", SPS_Ex - res.reconEx);
	fillAll2D("intermediateExIMvsSPS", SPS_Ex, res.intermediateEx);
	fillAll("MissingMass",res.missingmass);
	fillAll2D("intermediateExVSMissingMass", res.missingmass, res.intermediateEx);

	//intermediate CM
	fillAll("intermediatevcm_meas", res.intermediatevcm);
	//fillAll("intermediatevcm_expect", res.expected.vcm_intermediate);
	//fillAll("intermediatevcm_delta", res.intermediatevcm - res.expected.vcm_intermediate);
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
		//fillAll(f+"vcm_expect",expV);
		//fillAll(f+"vcm_delta",vcm-expV);
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
	fillAll("decay1_relLabAngle", res.relLabAngle_intfrag1);
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
	fillAll("decay2_relLabAngle", res.relLabAngle_frag2frag3);

	fillAll2D("decay2VSdecay1_relLabAngle", res.relLabAngle_intfrag1, res.relLabAngle_frag2frag3);

	fillAll2D("intermediatevcmVSfrag1vcm", res.frag1vcm, res.intermediatevcm);
	fillAll2D("intermediatekecmVSfrag1kecm", res.frag1kecm, res.intermediatekecm);
	fillAll2D("frag2vcmVSfrag3vcm", res.frag3vcm, res.frag2vcm);
	fillAll2D("frag2kecmVSfrag3kecm", res.frag3kecm, res.frag2kecm);
	fillAll2D("ecm1VSecm2", res.ecm2, res.ecm1);
	fillAll2D("ecm1deltaVSecm2delta", res.ecm2 - res.expected.Ecm2, res.ecm1 - res.expected.Ecm1);	
}

void MissMass_Mult2::FillGatedEventHistograms(double SPS_Ex){
	for(int i=0; i<6; i++) FillSelectGatedCaseHistograms(i, SPS_Ex);
}

void MissMass_Mult2::FillSelectGatedCaseHistograms(int caseNum, double SPS_Ex){

	if(caseNum < 0 || caseNum > 5){
		std::cerr << "caseNum out of range (" << caseNum << ")\n";
		return;
	}

	TString pn = permNames.at(caseNum);
	auto& res = caseResults[caseNum];

	if(caseResults[caseNum].permPasses){

		//lambdas to fill both specific permutation and "allCases"
		auto fillGated = [&](TString key, double x){
			groups_gated[pn]->Fill(key, x);
			groups_gated["allCases"]->Fill(key, x);
		};

		auto fillGated2D = [&](TString key, double x, double y){
			groups_gated[pn]->Fill(key, x, y);
			groups_gated["allCases"]->Fill(key, x, y);
		};

		//invariant mass and excitation energy histograms:
		fillGated("intermediateIM", res.intermediateIM);
		fillGated("intermediateEx", res.intermediateEx);
		fillGated("RecoilEx", res.reconEx);
		fillGated2D("RecoilEx_IMvsSPS", SPS_Ex, res.reconEx);
		fillGated("RecoilExDif", SPS_Ex - res.reconEx);
		fillGated2D("intermediateExIMvsSPS", SPS_Ex, res.intermediateEx);
		fillGated("MissingMass",res.missingmass);
		fillGated2D("intermediateExVSMissingMass", res.missingmass, res.intermediateEx);

		//intermediate CM
		fillGated("intermediatevcm_meas", res.intermediatevcm);
		//fillGated("intermediatevcm_expect", res.expected.vcm_intermediate);
		//fillGated("intermediatevcm_delta", res.intermediatevcm - res.expected.vcm_intermediate);
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
			//fillGated(f+"vcm_expect",expV);
			//fillGated(f+"vcm_delta",vcm-expV);
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
		//fillGated("ecm1_expect", res.expected.Ecm1);
		//fillGated("ecm1_delta", res.ecm1 - res.expected.Ecm1);
		fillGated2D("ecm1measVSintermediatethetacm", res.intermediatethetacm, res.ecm1);
		fillGated2D("ecm1measVSfrag1thetacm", res.frag1thetacm, res.ecm1);
		fillGated2D("ecm1measVSfrag2thetacm", res.frag2thetacm, res.ecm1);
		fillGated2D("ecm1measVSfrag3thetacm", res.frag3thetacm, res.ecm1);
		fillGated("decay1_VCM", std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1] + res.boost1[2]*res.boost1[2]));
		fillGated2D("decay1_VCM_TransverseVSLongitudinal", std::abs(res.boost1[2]), std::sqrt(res.boost1[0]*res.boost1[0] + res.boost1[1]*res.boost1[1]));
		fillGated("decay1_thetaCMsum", res.intermediatethetacm + res.frag1thetacm);
		fillGated("decay1_phiCMdiff", std::abs(res.intermediatephicm - res.frag1phicm));
		fillGated("decay1_relLabAngle", res.relLabAngle_intfrag1);
		fillGated("ecm2_meas", res.ecm2);
		//fillGated("ecm2_expect", res.expected.Ecm2);
		//fillGated("ecm2_delta", res.ecm2 - res.expected.Ecm2);
		fillGated2D("ecm2measVSintermediatethetacm", res.intermediatethetacm, res.ecm2);
		fillGated2D("ecm2measVSfrag1thetacm", res.frag1thetacm, res.ecm2);
		fillGated2D("ecm2measVSfrag2thetacm", res.frag2thetacm, res.ecm2);
		fillGated2D("ecm2measVSfrag3thetacm", res.frag3thetacm, res.ecm2);
		fillGated("decay2_VCM", std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1] + res.boost2[2]*res.boost2[2]));
		fillGated2D("decay2_VCM_TransverseVSLongitudinal", std::abs(res.boost2[2]), std::sqrt(res.boost2[0]*res.boost2[0] + res.boost2[1]*res.boost2[1]));
		fillGated("decay2_thetaCMsum", res.frag2thetacm + res.frag3thetacm);
		fillGated("decay2_phiCMdiff", std::abs(res.frag2phicm - res.frag3phicm));
		fillGated("decay2_relLabAngle", res.relLabAngle_frag2frag3);

		fillGated2D("decay2VSdecay1_relLabAngle", res.relLabAngle_intfrag1, res.relLabAngle_frag2frag3);

		fillGated2D("intermediatevcmVSfrag1vcm", res.frag1vcm, res.intermediatevcm);
		fillGated2D("intermediatekecmVSfrag1kecm", res.frag1kecm, res.intermediatekecm);
		fillGated2D("frag2vcmVSfrag3vcm", res.frag3vcm, res.frag2vcm);
		fillGated2D("frag2kecmVSfrag3kecm", res.frag3kecm, res.frag2kecm);
		fillGated2D("ecm1VSecm2", res.ecm2, res.ecm1);
		fillGated2D("ecm1deltaVSecm2delta", res.ecm2 - res.expected.Ecm2, res.ecm1 - res.expected.Ecm1);
	}
}

void MissMass_Mult2::FillPermCounter(bool gated){
	if(gated) hPermCounter_gated->Fill(CountPermPasses());
	else hPermCounter->Fill(CountPermPasses());
}

void MissMass_Mult2::FillSortedHisto(double SPS_Ex){

	ResultsMM sorted_caseResults[6];
	for(int i=0; i<6; i++) sorted_caseResults[i] = caseResults[i];

	std::sort(std::begin(sorted_caseResults), std::end(sorted_caseResults), [SPS_Ex](const ResultsMM& a, const ResultsMM& b) {
		double diffA = std::abs(SPS_Ex - a.reconEx);
		double diffB = std::abs(SPS_Ex - b.reconEx);
		return diffA < diffB;
	});

	hSortedIntermediateExIMvsSPS->Fill(SPS_Ex, sorted_caseResults[0].intermediateEx);
	hSortedIMRecEx->Fill(sorted_caseResults[0].reconEx);
	if(sorted_caseResults[0].intermediateEx >= -0.04 && sorted_caseResults[0].intermediateEx <= 0.04){
		hSortedIMRecEx_gate8Be->Fill(sorted_caseResults[0].reconEx);
	}
	if(sorted_caseResults[0].intermediateEx >= -0.8 && sorted_caseResults[0].intermediateEx <= 0.8){
		hSortedIMRecEx_gate5Li->Fill(sorted_caseResults[0].reconEx);
	}

	int permIndex = -1;
	TString perm = sorted_caseResults[0].permName;
	if(perm == "recon_f1a") permIndex = 0;
	else if(perm == "recon_f1b") permIndex = 1;
	else if(perm == "recon_f2a") permIndex = 2;
	else if(perm == "recon_f2b") permIndex = 3;
	else if(perm == "recon_f3a") permIndex = 4;
	else if(perm == "recon_f3b") permIndex = 5;

	hSortedPermutations->Fill(permIndex);
}

void MissMass_Mult2::FillSABREvsSPSHisto(double SPS_Ex, double SABREsumE){
	hSABRESumE_vs_ExSPS->Fill(SPS_Ex, SABREsumE);
}

void MissMass_Mult2::CloseAndWrite(){
	outfile->cd();
	hPermCounter->Write();
	hPermCounter_gated->Write();
	hSortedIntermediateExIMvsSPS->Write();
	hSortedPermutations->Write();
	hSortedIMRecEx->Write();
	hSortedIMRecEx_gate8Be->Write();
	hSortedIMRecEx_gate5Li->Write();
	hSABRESumE_vs_ExSPS->Write();
	if(outfile && outfile->IsOpen()){
		outfile->Write();
		outfile->Close();
	}
}

void MissMass_Mult2::SetExpectedCMValues(bool verbose){
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

int MissMass_Mult2::CountPermPasses(){
	int numPasses = 0;
	for(int i=0; i<6; i++) numPasses += caseResults[i].permPasses;
	return numPasses;
}

void MissMass_Mult2::ClearEventResults(){
	for(int i=0; i<6; i++) caseResults[i].Reset();
}
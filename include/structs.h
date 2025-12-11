#ifndef STRUCTS_H
#define STRUCTS_H

static const int numwedges = 8;
static const int numrings = 16;

static const float DEGRAD=0.017453293;

static const std::pair<int,int> offsets[] = {
	{112,40},	//detector0 {ringOffset,wedgeOffset}
	{96,32},	//detector1 {ringOffset,wedgeOffset}
	{80,16},	//detector2 {ringOffset,wedgeOffset}
	{64,24},	//detector3 {ringOffset,wedgeOffset}
	{48,0}		//detector4 {ringOffset,wedgeOffset}
};

struct PHYSDATA { double e, theta, phi; };
struct SABREDATA { int detectorIndex=-666, particleIndex=-666; double theta, phi, ringEnergy, wedgeEnergy, localx, localy; int ring, wedge; };

struct Nucleus {
	int A;
	TString sym;
	double massMeV;

	void SetA(int newA){
		A = newA;
	}

	void SetSym(TString newsym){
		sym = newsym;
	}

	void SetMassMeV(double newmassMeV){
		massMeV = newmassMeV;
	}

	void SetAll(int newA, TString newsym, double newmassMeV){
		SetA(newA);
		SetSym(newsym);
		SetMassMeV(newmassMeV);
	}
};

struct Reaction2{
	Nucleus beam;
	Nucleus target;
	Nucleus ejectile;
	Nucleus recoil;

	double beamEnergy=0.;
	double recoilExE=0.;

	TString ToString(){
		TString retval = Form("%d%s(%d%s,%d%s)%d%s at E=%f to ExE=%f",target.A,target.sym.Data(),beam.A,beam.sym.Data(),ejectile.A,ejectile.sym.Data(),recoil.A,recoil.sym.Data(),beamEnergy,recoilExE);
		return retval;
	};
};

struct Reaction3{
	Nucleus beam;
	Nucleus target;
	Nucleus ejectile;
	Nucleus recoil;
	Nucleus breakup1;
	Nucleus breakup2;

	double beamEnergy=0.;
	double recoilExE=0.;

	TString ToString(){
		TString retval = Form("%d%s(%d%s,%d%s)%d%s at E=%f to ExE=%f",target.A,target.sym.Data(),beam.A,beam.sym.Data(),ejectile.A,ejectile.sym.Data(),recoil.A,recoil.sym.Data(),beamEnergy,recoilExE);
		return retval;
	};
};

#endif//STRUCTS_H
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

int ReadSRIMTransmitIntoTTree(TString infilename, TString outfilename){

	const int numheaderlines = 12;

	std::ifstream fin(infilename.Data());
	if(!fin.is_open()){
		std::cerr << "Error: Cannot open " << infilename << std::endl;
		return -1;
	}

	TFile *outputfile = new TFile(outfilename, "RECREATE");

	TString treetitle = outfilename;
	int lastslash = treetitle.Last('/');
	if(lastslash>=0) treetitle = treetitle(lastslash+1, treetitle.Length()-(lastslash+1));
	treetitle.ReplaceAll(".root","");


	TTree *tree = new TTree("SRIMTransmit", treetitle);

	char type;
	int ionNum, atomNum;
	double energy, x, y, z, cosx, cosy, cosz;

	//branches:
	tree->Branch("type", &type, "type/C");
	tree->Branch("ionNum", &ionNum, "ionNum/I");
	tree->Branch("atomNum", &atomNum, "atomNum/I");
	tree->Branch("energy", &energy, "energy/D");
	tree->Branch("x", &x, "x/D");
	tree->Branch("y", &y, "y/D");
	tree->Branch("z", &z, "z/D");
	tree->Branch("cosx", &cosx, "cosx/D");
	tree->Branch("cosy", &cosy, "cosy/D");
	tree->Branch("cosz", &cosz, "cosz/D");

	std::string line;


	int entry = 0;
	while(std::getline(fin,line)){
		std::istringstream iss(line);
		if(iss >> type >> ionNum >> atomNum >> energy >> x >> y >> z >> cosx >> cosy >> cosz){
			entry += 1;
			//if(ionNum != entry) std::cout << "filling entry number " << entry << " which is ionNum = " << ionNum << "\n";
			tree->Fill();
			if(ionNum%1000==0) std::cout << "ionnum = " << ionNum << ", tree entry = " << tree->GetEntries() << std::endl;

		}
	}

	//tree->Write();
	std::cout << "TTree with name SRIMTransmit and title " << treetitle << " and " << tree->GetEntries() << " entries written to " << outfilename << "\n";

	outputfile->Write();
	outputfile->Close();
	fin.close();
	return 0;
}
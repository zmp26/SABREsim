void CompareFiles(const char* SABRE_file, const char* kin4mc_file){
	TFile *sabrefile = TFile::Open(SABRE_file, "READ");
	TFile *kin4mcfile = TFile::Open(kin4mc_file, "READ");

	TTree *sabretree = (TTree*)sabrefile->Get("mult3");
	TTree *kin4mctree = (TTree*)kin4mcfile->Get("mult3");

	sabretree->AddFriend(kin4mctree, "kin");

	sabretree->Draw("frag1thetacm:kin.frag1thetacm");
}
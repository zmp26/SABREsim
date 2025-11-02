#include "RootWriter.h"
#include <cstring>

RootWriter::RootWriter(const std::string& filename){
	file_ = new TFile(filename.c_str(),"RECREATE");
	tree_ = new TTree("SABREsim", "SABREsim");

	tree_->Branch("kin_e", kin_e_, "kin_e[4]/D");
	tree_->Branch("kin_theta", kin_theta_, "kin_theta[4]/D");
	tree_->Branch("kin_phi", kin_phi_, "kin_phi[4]/D");

	tree_->Branch("reactionOrigin",reactionOrigin_, "reactionOrigin[3]/D");

	tree_->Branch("numHits", &numHits_, "numHits/I");
	tree_->Branch("particleID", particleID_, "particleID[4]/I");
	tree_->Branch("detectorID", detectorID_, "detectorID[4]/I");
	tree_->Branch("ringChannel", ringChannel_, "ringChannel[4]/I");
	tree_->Branch("wedgeChannel", wedgeChannel_, "wedgeChannel[4]/I");
	tree_->Branch("localRing", localRing_, "localRing[4]/I");
	tree_->Branch("localWedge", localWedge_, "localWedge[4]/I");
	tree_->Branch("ringEnergy", ringEnergy_, "ringEnergy[4]/D");
	tree_->Branch("wedgeEnergy", wedgeEnergy_, "wedgeEnergy[4]/D");
	tree_->Branch("localx", localx_, "localx[4]/D");
	tree_->Branch("localy", localy_, "localy[4]/D");

	ResetEvent();

	//default metadata values:
	detmcVersion_ = 0;
	reaction_ = "";
	beamEnergyMeV_ = -666.;
	beamSpotProfile_ = "";
	beamSpotParX_ = 0.;
	beamSpotParY_ = 0.;
}

RootWriter::~RootWriter() {
	WriteAndClose();
}

void RootWriter::ResetEvent() {
	std::fill_n(kin_e_, 4, -666.);
	std::fill_n(kin_theta_, 4, -666.);
	std::fill_n(kin_phi_, 4, -666.);

	std::fill_n(reactionOrigin_, 3, -666.);

	numHits_ = 0;
	std::fill_n(particleID_, 4, -666);
	std::fill_n(detectorID_, 4, -666);
	std::fill_n(ringChannel_, 4, -666);
	std::fill_n(wedgeChannel_, 4, -666);
	std::fill_n(localRing_, 4, -666);
	std::fill_n(localWedge_, 4, -666);
	std::fill_n(ringEnergy_, 4, -666.);
	std::fill_n(wedgeEnergy_, 4, -666.);
	std::fill_n(localx_, 4, -666.);
	std::fill_n(localy_, 4, -666.);

}

void RootWriter::SetKinematics(int index, double e, double theta, double phi){
	if(index < 0 || index >=4 ) return;
	kin_e_[index] = e;
	kin_theta_[index] = theta;
	kin_phi_[index] = phi;
}

void RootWriter::SetReactionOrigin(double x, double y, double z){
	reactionOrigin_[0] = x;
	reactionOrigin_[1] = y;
	reactionOrigin_[2] = z;
}

void RootWriter::AddHit(int hitIndex,
						int particleID,
						int detectorID,
						int ringChannel,
						int wedgeChannel,
						int localRing,
						int localWedge,
						double ringEnergy,
						double wedgeEnergy,
						double localx,
						double localy)
{
	if(hitIndex < 0 || hitIndex >= 4) return;

	particleID_[hitIndex] = particleID;
	detectorID_[hitIndex] = detectorID;
	ringChannel_[hitIndex] = ringChannel;
	wedgeChannel_[hitIndex] = wedgeChannel;
	localRing_[hitIndex] = localRing;
	localWedge_[hitIndex] = localWedge;
	ringEnergy_[hitIndex] = ringEnergy;
	wedgeEnergy_[hitIndex] = wedgeEnergy;
	localx_[hitIndex] = localx;
	localy_[hitIndex] = localy;
}

void RootWriter::FillEvent(){
	for(const auto& e : ringEnergy_){
		if(e > 0) numHits_++;
	}
	tree_->Fill();
	ResetEvent();
}

//metadata:
void RootWriter::WriteMetaData(){
	if(!file_ || !file_->IsOpen()) return;

	file_->cd();

	TParameter<int>("detmcVersion",detmcVersion_).Write();
	TParameter<std::string>("reaction",reaction_).Write();
	TParameter<double>("beamEnergyMeV",beamEnergyMeV_).Write();
	TParameter<std::string>("beamSpotProfile",beamSpotProfile_).Write();
	TParameter<double>("beamSpotParX",beamSpotParX_).Write();
	TParameter<double>("beamSpotParY",beamSpotParY_).Write();
}

void RootWriter::WriteAndClose(){
	if(file_ && file_->IsOpen()){
		file_->cd();
		WriteMetaData();
		tree_->Write();
		file_->Close();
		delete file_;
		file_ = nullptr;
	}
}
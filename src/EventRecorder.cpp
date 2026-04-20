#include "EventRecorder.h"
#include "TNamed.h"
#include "TParameter.h"
#include <algorithm>

EventRecorder::EventRecorder(const std::string& filename)
	: file_(nullptr), tree_(nullptr), numHits_(0),
	inputfile_(""), detmcVersion_(0), reaction_(""), beamEnergyMeV_(-666.),
	beamSpotProfile_(""), beamSpotParX_(0.), beamSpotParY_(0.)
{
	file_ = new TFile(filename.c_str(), "RECREATE");
	tree_ = new TTree("SABREsim", "SABREsim");

	tree_->Branch("kin_e", kin_e_, "kin_e[4]/D");
	tree_->Branch("kin_theta", kin_theta_, "kin_theta[4]/D");
	tree_->Branch("kin_phi", kin_phi_, "kin_phi[4]/D");
	tree_->Branch("reactionOrigin", reactionOrigin_, "reactionOrigin[3]/D");

	tree_->Branch("EjInSPS", &EjInSPS, "EjInSPS/O");
	tree_->Branch("SPSEnergy", &SPSEnergy, "SPSEnergy/D");
	tree_->Branch("SPSTheta", &SPSTheta, "SPSTheta/D");
	tree_->Branch("SPSPhi", &SPSPhi, "SPSPhi/D");

	tree_->Branch("numHits", &numHits_, "numHits/I");
	tree_->Branch("particleID", particleID_, "particleID[4]/I");
	tree_->Branch("detectorID", detectorID_, "detectorID[4]/I");
	tree_->Branch("ringChannel", ringChannel_, "ringChannel[4]/I");
	tree_->Branch("wedgeChannel", wedgeChannel_, "wedgeChannel[4]/I");
	tree_->Branch("localRing", localRing_, "localRing[4]/I");
	tree_->Branch("localWedge", localWedge_, "localWedge[4]/I");
	tree_->Branch("ringEnergy", ringEnergy_, "ringEnergy[4]/D");
	tree_->Branch("wedgeEnergy", wedgeEnergy_, "wedgeEnergy[4]/D");
	tree_->Branch("ringTheta", ringTheta_, "ringTheta[4]/D");
	tree_->Branch("wedgePhi", wedgePhi_, "wedgePhi[4]/D");
	tree_->Branch("localx", localx_, "localx[4]/D");
	tree_->Branch("localy", localy_, "localy[4]/D");

	ResetEvent();
}

EventRecorder::~EventRecorder(){
	WriteAndClose();
}

void EventRecorder::ResetEvent(){
	numHits_ = 0;
	hits_.clear();

	std::fill_n(kin_e_, 4, -666.);
	std::fill_n(kin_theta_, 4, -666.);
	std::fill_n(kin_phi_, 4, -666.);
	std::fill_n(reactionOrigin_, 3, -666.);

	EjInSPS = false;
	SPSEnergy = -666.;
	SPSTheta = -666.;
	SPSPhi = -666.;

	std::fill_n(particleID_, 4, -666);
	std::fill_n(detectorID_, 4, -666);
	std::fill_n(ringChannel_, 4, -666);
	std::fill_n(wedgeChannel_, 4, -666);
	std::fill_n(localRing_, 4, -666);
	std::fill_n(localWedge_, 4, -666);
	std::fill_n(ringEnergy_, 4, -666.);
	std::fill_n(wedgeEnergy_, 4, -666.);
	std::fill_n(ringTheta_, 4, -666.);
	std::fill_n(wedgePhi_, 4, -666.);
	std::fill_n(localx_, 4, -666.);
	std::fill_n(localy_, 4, -666.);
}

// void EventRecorder::SetNumParticles(int n){
// 	nParticles_ = std::min(std::max(n,1),4);
// }

void EventRecorder::SetKinematics(int index, double e, double theta, double phi){
	if(index < 0 || index >= 4) return;
	kin_e_[index] = e;
	kin_theta_[index] = theta;
	kin_phi_[index] = phi;
}

void EventRecorder::SetReactionOrigin(double x, double y, double z){
	reactionOrigin_[0] = x;
	reactionOrigin_[1] = y;
	reactionOrigin_[2] = z;
}

void EventRecorder::UpdateSPS(bool EjectileInSPS, double SPSE, double SPSTh, double SPSPh){
	EjInSPS = EjectileInSPS;
	SPSEnergy = SPSE;
	SPSTheta = SPSTh;
	SPSPhi = SPSPh;
}

void EventRecorder::AddHit(const Hit& hit){
	hits_.push_back(hit);
}

void EventRecorder::FillEvent(){


	int n = std::min((int)hits_.size(), 4);
	for(int i=0; i<n; i++){
		const Hit& h = hits_[i];
		particleID_[i] = h.particleID;
		detectorID_[i] = h.detectorID;
		ringChannel_[i] = h.ringChannel;
		wedgeChannel_[i] = h.wedgeChannel;
		localRing_[i] = h.localRing;
		localWedge_[i] = h.localWedge;
		ringEnergy_[i] = h.ringEnergy;
		wedgeEnergy_[i] = h.wedgeEnergy;
		ringTheta_[i] = h.ringTheta;
		wedgePhi_[i] = h.wedgePhi;
		localx_[i] = h.localx;
		localy_[i] = h.localy;
		if(h.ringEnergy > 0) numHits_ += 1;
	}

	tree_->Fill();
	ResetEvent();
}

void EventRecorder::WriteMetaData(){
	if(!file_ || !file_->IsOpen()) return;

	file_->cd();
	TNamed("inputFile", inputfile_.c_str()).Write();
	TParameter<int>("detmcVersion",detmcVersion_).Write();
	TNamed("reaction", reaction_.c_str()).Write();
	TParameter<double>("beamEnergyMeV",beamEnergyMeV_).Write();
	TNamed("beamSpotProfile",beamSpotProfile_.c_str()).Write();
	TParameter<double>("beamSpotParX", beamSpotParX_).Write();
	TParameter<double>("beamSpotParY", beamSpotParY_).Write();
}

void EventRecorder::WriteAndClose(){
	if(file_ && file_->IsOpen()){
		file_->cd();
		WriteMetaData();
		tree_->Write("",TObject::kOverwrite);
		file_->Close();
		delete file_;
		file_ = nullptr;
	}
}
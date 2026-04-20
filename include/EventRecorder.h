#ifndef EVENTRECORDER_H
#define EVENTRECORDER_H

#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "structs.h"

class EventRecorder{
public:
	EventRecorder(const std::string& filename);
	~EventRecorder();

	void ResetEvent();
	void SetKinematics(int index, double e, double theta, double phi);
	void SetReactionOrigin(double x, double y, double z);
	void UpdateSPS(bool EjInSPS, double SPSE, double SPSTh, double SPSPh);
	void AddHit(const Hit& hit);
	void FillEvent();

	void WriteMetaData();
	void WriteAndClose();

	void SetInputFile(const std::string& f) { inputfile_ = f; }
	void SetDetMCVersion(int v) { detmcVersion_ = v; }
	void SetReaction(const std::string& r) { reaction_ = r; }
	void SetBeamEnergyMeV(double e) { beamEnergyMeV_ = e; }
	void SetBeamSpotProfile(const std::string& s) { beamSpotProfile_ = s; }
	void SetBeamSpotParameters(double x, double y) { beamSpotParX_ = x; beamSpotParY_ = y; }
	void SetBeamSpotParX(double x) { beamSpotParX_ = x; }
	void SetBeamSpotParY(double y) { beamSpotParY_ = y; }

private:
	TFile *file_;
	TTree *tree_;

	double kin_e_[4];
	double kin_theta_[4];
	double kin_phi_[4];
	double reactionOrigin_[3];

	bool EjInSPS;
	double SPSEnergy, SPSTheta, SPSPhi;

	int numHits_;
	int particleID_[4];
	int detectorID_[4];
	int ringChannel_[4];
	int wedgeChannel_[4];
	int localRing_[4];
	int localWedge_[4];
	double ringEnergy_[4];
	double wedgeEnergy_[4];
	double ringTheta_[4];
	double wedgePhi_[4];
	double localx_[4];
	double localy_[4];

	std::vector<Hit> hits_;

	std::string inputfile_;
	int detmcVersion_;
	std::string reaction_;
	double beamEnergyMeV_;
	std::string beamSpotProfile_;
	double beamSpotParX_;
	double beamSpotParY_;

};

#endif//EVENTRECORDER_H
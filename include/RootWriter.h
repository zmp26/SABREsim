#ifndef ROOTWRITER_H
#define ROOTWRITER_H


#include <vector>
#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"


class RootWriter {
public:
	RootWriter(const std::string& filename);
	~RootWriter();

	void ResetEvent();

	//kinematics
	void SetKinematics(int index, double e, double theta, double phi);

	//reaction origin
	void SetReactionOrigin(double x, double y, double z);

	//sabre detection info
	void AddHit(int hitIndex,
				int particleID,
				int detectorID,
				int ringChannel,
				int wedgeChannel,
				int localRing,
				int localWedge,
				double ringEnergy,
				double wedgeEnergy,
				double localx,
				double localy);

	void FillEvent();
	void WriteAndClose();

	//metadata:
	void Set_detmc(int version); //2,3,4
	void SetReaction(const std::string& reaction);
	void SetBeamEnergyMeV(double energy);
	void SetBeamSpotProfile(const std::string& profile);
	void SetBeamSpotParameters(double x, double y);

private:
	TFile *file_;
	TTree *tree_;

	//kinematics
	double kin_e_[4];
	double kin_theta_[4];
	double kin_phi_[4];

	//reaction origin
	double reactionOrigin_[3];

	//int sabre hits
	int numHits_;
	int particleID_[4];
	int detectorID_[4];
	int ringChannel_[4];
	int wedgeChannel_[4];
	int localRing_[4];
	int localWedge_[4];
	double ringEnergy_[4];
	double wedgeEnergy_[4];
	double localx_[4];
	double localy_[4];

	//metadata:
	int detmcVersion_;
	std::string reaction_;
	double beamEnergyMeV_;
	std::string beamSpotProfile_;
	double beamSpotParX_;
	double beamSpotParY_;

	void WriteMetaData();
};


#endif//ROOTWRITER_H
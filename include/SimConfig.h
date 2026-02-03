#ifndef SIMCONFIG_H
#define SIMCONFIG_H

#include <string>
#include <vector>
#include <unordered_map>
#include "MassTable.h"
#include "TString.h"

struct NucleusConfig {
	int A = 0;
	TString symbol = "";
	double mass = 0.0;
};

class SimConfig{
public:
	explicit SimConfig(const std::string& filename);
	~SimConfig();
	bool Parse();

	int GetDetMCVersion() const { return detmc_version_; }
	const std::string& GetInFile() { return infile_; }
	const std::string& GetDetFile() { return detfile_; }
	const std::string& GetTreeFile() { return treefile_; }
	const std::string& GetHistoFile() { return histofile_; }

	const std::string& GetTargetLoss(int i) const { return targetLoss_par_.at(i-1); }
	const std::string& GetDeadLayerLoss(int i) const { return deadLayerLoss_par_.at(i-1); }
	const std::string& GetTargetStraggler(int i) const { return targetStraggle_par_.at(i-1); }

	const std::string& GetBeamProfile() const { return beam_profile_; }
	double GetBeamParX() const { return beam_parX_; }
	double GetBeamParY() const { return beam_parY_; }
	double GetBeamOffsetX() const { return beam_offsetX_; }
	double GetBeamOffsetY() const { return beam_offsetY_; }

	const std::string& GetReaction() const { return reaction_; }
	double GetBeamEnergy() const { return beam_energy_; }
	double GetRecoilExcitationEnergy() const { return recoil_excitation_energy_; }

	bool GetStraggleEnabled(int i) const { return enableStraggle_par_.at(i-1); }
	std::string GetStraggle(int i) const { return targetStraggle_par_.at(i-1); }

	const NucleusConfig& GetBeam() const { return beam_; }
	const NucleusConfig& GetTarget() const { return target_; }
	const NucleusConfig& GetEjectile() const { return ejectile_; }
	const NucleusConfig& GetRecoil() const { return recoil_; }

	const std::vector<NucleusConfig>& GetBreakups() const { return breakups_; }

private:
	MassTable *masstable;

	std::string filename_;

	int detmc_version_;
	std::string infile_, detfile_, treefile_, histofile_;

	std::vector<std::string> targetLoss_par_;

	//double straggleMu_, straggleSigma_, straggleLambda_;
	std::vector<bool> enableStraggle_par_;
	std::vector<std::string> targetStraggle_par_;

	std::vector<std::string> deadLayerLoss_par_;

	std::string beam_profile_;
	double beam_parX_, beam_parY_;
	double beam_offsetX_, beam_offsetY_;

	std::string reaction_;
	double beam_energy_;
	double recoil_excitation_energy_;

	std::string Trim(const std::string& s);

	NucleusConfig beam_;
	NucleusConfig target_;
	NucleusConfig ejectile_;
	NucleusConfig recoil_;
	std::vector<NucleusConfig> breakups_;
};


#endif//SIMCONFIG_H
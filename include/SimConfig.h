#ifndef SIMCONFIG_H
#define SIMCONFIG_H

#include <string>
#include <vector>
#include <unordered_map>

class SimConfig{
public:
	explicit SimConfig(const std::string& filename);
	bool Parse();

	int GetDetMCVersion() const { return detmc_version_; }
	const std::string& GetInFile() { return infile_; }
	const std::string& GetDetFile() { return detfile_; }
	const std::string& GetTreeFile() { return treefile_; }
	const std::string& GetHistoFile() { return histofile_; }

	const std::string& GetTargetLoss(int i) const { return targetLoss_par_.at(i-1); }
	const std::string& GetDeadLayerLoss(int i) const { return deadLayerLoss_par_.at(i-1); }

	const std::string& GetBeamProfile() const { return beam_profile_; }
	double GetBeamParX() const { return beam_parX_; }
	double GetBeamParY() const { return beam_parY_; }
	double GetBeamOffsetX() const { return beam_offsetX_; }
	double GetBeamOffsetY() const { return beam_offsetY_; }

	const std::string& GetReaction() const { return reaction_; }
	double GetBeamEnergy() const { return beam_energy_; }
	double GetRecoilExcitationEnergy() const { return recoil_excitation_energy_; }

	double GetStraggleMu() const { return straggleMu_; }
	double GetStraggleSigma() const { return straggleSigma_; }
	double GetStraggleLambda() const { return straggleLambda_; }


private:
	std::string filename_;

	int detmc_version_;
	std::string infile_, detfile_, treefile_, histofile_;

	std::vector<std::string> targetLoss_par_;
	double straggleMu_, straggleSigma_, straggleLambda_;
	std::vector<std::string> deadLayerLoss_par_;

	std::string beam_profile_;
	double beam_parX_, beam_parY_;
	double beam_offsetX_, beam_offsetY_;

	std::string reaction_;
	double beam_energy_;
	double recoil_excitation_energy_;

	std::string Trim(const std::string& s);
};


#endif//SIMCONFIG_H
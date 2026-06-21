#ifndef HISTOMANAGER_H
#define HISTOMANAGER_H

#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"
#include <map>
#include <string>
#include <memory>
#include <vector>

class HistoManager {
private:
	std::map<std::string, std::unique_ptr<TH1>> hMap;
	TDirectory* dir;

public:
	HistoManager(TDirectory* targetDir) : dir(targetDir) {}

	void Register1D(std::string key, std::string title, int bins, double min, double max);
	void Register2D(std::string key, std::string title, int xbins, double xmin, double xmax, double ybins, double ymin, double ymax);

	//void RegisterNodeHistograms(std::string prefix, std::string nodeName);
	void RegisterHistograms(const KinematicNode* node, std::string perm);

	void Fill(std::string key, double x);
	void Fill(std::string key, double x, double y);
};

#endif//HISTOMANAGER_H
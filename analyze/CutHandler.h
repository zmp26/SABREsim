#ifndef CUT_HANDLER_H
#define CUT_HANDLER_H

#include "TFile.h"
#include "TCutG.h"
#include "TString.h"
#include <vector>
#include <map>
#include <iostream>

struct Cut1D {
	TString name;
	double minVal;
	double maxVal;
	bool isThreshold;//true only if checking above single value
	bool isUpperLimit;//true only if checking below single value
};

class CutHandler {
private:
	std::vector<Cut1D> fCuts1D;
	std::vector<TCutG*> fCuts2D;

public:
	CutHandler() = default;
	~CutHandler();

	//1D cut registration
	void AddCut1D(const TString& name, double min, double max, bool isThreshold=false, bool isUpperLimit=false);

	//2D cut registration
	bool LoadCut2D(const TString& fileName, const TString& cutName);
	bool LoadAllCuts2D(const TString& fileName, const std::vector<TString>& cutNames);
	void AddDirectCut2D(TCutG* cut);

	bool CheckCut1D(const TString& name, double value) const;
	bool CheckCut2D(const TString& name, double x, double y) const;

	//batch check all 1D/2D gates
	bool PassAll1D(const std::map<TString, double>& values) const;
	bool PassAll2D(const std::map<TString, std::pair<double,double>>& eventPoints) const;

	void PrintCuts() const;
	const std::vector<TCutG*>& GetCuts2D() const { return fCuts2D; }
	const std::vector<Cut1D>& GetCuts1D() const { return fCuts1D; }
};

#endif//CUT_HANDLER_H
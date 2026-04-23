#ifndef PERMHISTO_H
#define PERMHISTO_H

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TDirectory.h"
#include <map>

class permHisto {
public:
	permHisto(TString permName, TDirectory* targetDir);

	~permHisto(){};

	void Fill(TString key, double x);
	void Fill(TString key, double x, double y);

	TH1* Get(TString key) { return hMap[key]; }

private:
	std::map<TString, TH1*> hMap;

	void Register1D(TString permName, TString key, TString title, int xbins, double xmin, double xmax);
	void Register2D(TString permName, TString key, TString title, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax);
};

#endif//PERMHISTO_H
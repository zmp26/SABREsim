#ifndef PERMHISTOCORRELATIONMULT3_H
#define PERMHISTOCORRELATIONMULT3_H

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TDirectory.h"
#include <map>
#include <vector>
#include <string>

class permHistoCorrelation_mult3 {
public:
    permHistoCorrelation_mult3(const std::vector<TString>& permNames, TDirectory* targetDir);
    ~permHistoCorrelation_mult3() {};

    void FillCorner(const std::vector<double>& chi2Values);

private:
    std::map<TString, TH2D*> h2Map;
    std::map<TString, TH1D*> h1Map;
    std::vector<TString> fPermNames;
};

#endif //PERMHISTOCORRELATIONMULT3_H
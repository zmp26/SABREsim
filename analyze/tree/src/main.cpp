#include <iostream>
#include <vector>
#include <algorithm>
#include "EventMapper.h"
#include "HistoManager.h"
#include "TFile.h"

std::string GetPermName(const std::vector<int>& p){
    std::string name = "";
    for(int idx : p) name += std::to_string(idx);
    return name;
}

int main(){

    TFile *outfile = new TFile("filenamehere.root","RECREATE");
    HistoManager histoManager(&outfile);
    EventMapper mapper;

    if(!mapper.Initialize("configfilehere.txt")) return -1;

    std::vector<int> perms = {0, 1, 2};

    int permIndex = 0;
    do{

        std::string permName = GetPermName(p);
        histoManager.RegisterHistograms(mapper.GetRoot(), permName);

    } while(std::next_permutation(p.begin(), p.end()));

    //main event loop below:
    while(){

    }

}
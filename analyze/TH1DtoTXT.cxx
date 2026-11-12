#include <iostream>
#include <fstream>
#include "TH1D.h"

void TH1DtoTXT(TH1D* h, const char* filename = "output.txt") {
    if (!h) {
        std::cerr << "Error: Histogram pointer is null!" << std::endl;
        return;
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write header
    //outfile << "# Bin\tBinCenter\tBinContent\tBinError\n";

    int nbins = h->GetNbinsX();
    
    // Loop including underflow (0) and overflow (nbins + 1)
    for (int i = 0; i <= nbins + 1; ++i) {
        // outfile << i << "\t"
        //         << h->GetBinCenter(i) << "\t"
        //         << h->GetBinContent(i) << "\t"
        //         << h->GetBinError(i) << "\n";
        outfile << h->GetBinCenter(i) << ", " << h->GetBinContent(i) << "\n";
    }

    outfile.close();
    std::cout << "Successfully exported " << h->GetName() << " to " << filename << std::endl;
}

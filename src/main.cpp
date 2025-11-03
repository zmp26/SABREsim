using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include "ConsoleColorizer.h"
#include "SABREsim.h"
#include <TString.h>

int main(int argc, char * argv[]){

	if(argc == 2 && (std::string(argv[1]) == "help" || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")){
		// ConsoleColorizer::PrintWhite("Usage:\n\t");
		// ConsoleColorizer::PrintYellow("bin/SABREsim X path/to/kinXmcFile.out path/to/detOutputFile.det\n\n");
		// ConsoleColorizer::PrintWhite("Where ");
		// ConsoleColorizer::PrintYellow("X = 2,3,4 ");
		// ConsoleColorizer::PrintWhite("for kin2mc, kin3mc, kin4mc input  files\nAnd ");
		// ConsoleColorizer::PrintYellow("kinInputFile.out ");
		// ConsoleColorizer::PrintWhite("is the filename (or path to) the kinXmc output file\nAnd ");
		// ConsoleColorizer::PrintYellow("detOutputFile.det ");
		// ConsoleColorizer::PrintWhite("is the filename (or path to) the detection output file (what this code writes to)\n\nExample command:\n");
		// ConsoleColorizer::PrintYellow("\tbin/SABREsim 2 kinmc/TEST.out det/TEST.det\n\n");
		// return 1;

		ConsoleColorizer::PrintWhite("Usage:\n\t");
		ConsoleColorizer::PrintYellow("bin/SABREsim /path/to/config/file.conf\n\n");
		ConsoleColorizer::PrintWhite("Where ");
		ConsoleColorizer::PrintYellow("/path/to/config/file.conf ");
		ConsoleColorizer::PrintWhite("is the filename (or path to) the appropriate config file\n\n");
		return 1;
	}

	if(argc!=2){
		ConsoleColorizer::PrintRed("\nError! Please provide input as command line arguments!\nExpected command of form: ");
		ConsoleColorizer::PrintYellow("bin/SABREsim /path/to/config/file.conf\n");
		ConsoleColorizer::PrintRed("Where ");
		ConsoleColorizer::PrintYellow("/path/to/config/file.conf ");
		ConsoleColorizer::PrintRed("is the filename (or path to) the appropriate config file\n\n");
		return 1;
	}


	SABREsim sim(argv[1]);


	std::cout << "\n\n";
	ConsoleColorizer::PrintPurple("███████╗ █████╗ ██████╗ ██████╗ ███████╗███████╗██╗███╗   ███╗\n");
	ConsoleColorizer::PrintPurple("██╔════╝██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔════╝██║████╗ ████║\n");
	ConsoleColorizer::PrintPurple("███████╗███████║██████╔╝██████╔╝█████╗  ███████╗██║██╔████╔██║\n");
	ConsoleColorizer::PrintPurple("╚════██║██╔══██║██╔══██╗██╔══██╗██╔══╝  ╚════██║██║██║╚██╔╝██║\n");
	ConsoleColorizer::PrintPurple("███████║██║  ██║██████╔╝██║  ██║███████╗███████║██║██║ ╚═╝ ██║\n");
	ConsoleColorizer::PrintPurple("╚══════╝╚═╝  ╚═╝╚═════╝ ╚═╝  ╚═╝╚══════╝╚══════╝╚═╝╚═╝     ╚═╝\n");
                                                              

	std::cout << std::endl;
	TString msg = Form("Kin%dmc selected!\n\n",sim.GetKinX());
	ConsoleColorizer::PrintGreen(msg.Data());

	msg = Form("phys file:\t%s\n\n", sim.GetKinInputFilename().c_str());
	ConsoleColorizer::PrintGreen(msg.Data());

	msg = Form("det file:\t%s\n", sim.GetDetOutputFilename().c_str());
	ConsoleColorizer::PrintGreen(msg.Data());

	msg = Form("tree file:\t%s\n", sim.GetTreeFilename().c_str());
	ConsoleColorizer::PrintGreen(msg.Data());

	msg = Form("histo file:\t%s\n\n", sim.GetHistoFilename().c_str());
	ConsoleColorizer::PrintGreen(msg.Data());

	try{
		if(!sim.GetFailState()){
			sim.Run();
			return 0;
		} else {
			ConsoleColorizer::PrintRed("SABREsim fail state set to true! Double check kinX selection in config file and make sure it is 2, 3, or 4!\n\n");
			return 1;
		}

	} catch(const std::exception& e){
		ConsoleColorizer::PrintRed(Form("Exception: %s\n",e.what()));
		return 1;
	}
	catch(...){
		ConsoleColorizer::PrintRed("Unknown exception occured during simulation!\n");
		return 1;
	}
}
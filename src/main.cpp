using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include "ConsoleColorizer.h"
#include "SABREsim.h"

static const std::pair<int, int> offsets[] = {
	{112,40},	//detector0
	{96,32},	//detector1
	{80,16},	//detector2
	{64,24},	//detector3
	{48,0}		//detector4
};

int main(int argc, char * argv[]){

	if(argc == 2 && (std::string(argv[1]) == "help" || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")){
		ConsoleColorizer::PrintWhite("Usage:\n\t");
		ConsoleColorizer::PrintYellow("bin/SABREsim X path/to/kinXmcFile.out path/to/detOutputFile.det\n\n");
		ConsoleColorizer::PrintWhite("Where ");
		ConsoleColorizer::PrintYellow("X = 2,3,4 ");
		ConsoleColorizer::PrintWhite("for kin2mc, kin3mc, kin4mc input  files\nAnd ");
		ConsoleColorizer::PrintYellow("kinInputFile.out ");
		ConsoleColorizer::PrintWhite("is the filename (or path to) the kinXmc output file\nAnd ");
		ConsoleColorizer::PrintYellow("detOutputFile.det ");
		ConsoleColorizer::PrintWhite("is the filename (or path to) the detection output file (what this code writes to)\n\nExample command:\n");
		ConsoleColorizer::PrintYellow("\tbin/SABREsim 2 kinmc/TEST.out det/TEST.det\n\n");
		return 1;
	}

	if(argc!=4){
		ConsoleColorizer::PrintRed("\nError! Please provide input as command line arguments!\nExpected arguments in order of: ");
		ConsoleColorizer::PrintYellow("X kinXmcFile.out detOutputFile.det\n");
		ConsoleColorizer::PrintRed("Where ");
		ConsoleColorizer::PrintYellow("X = 2,3,4 ");
		ConsoleColorizer::PrintRed("for kin2mc, kin3mc, kin4mc input  files\nAnd ");
		ConsoleColorizer::PrintYellow("kinInputFile.out ");
		ConsoleColorizer::PrintRed("is the filename (or path to) the kinXmc output file\nAnd ");
		ConsoleColorizer::PrintYellow("detOutputFile.det ");
		ConsoleColorizer::PrintRed("is the filename (or path to) the detection output file (what this code writes to)\n\nExample command:\n");
		ConsoleColorizer::PrintYellow("\tbin/SABREsim 2 kinmc/TEST.out det/TEST.det\n\n");
		return 1;
	}

	int kinX = std::stoi(argv[1]);
	if(kinX != 2 && kinX !=3 && kinX !=4){
		ConsoleColorizer::PrintRed("Error: Invalid kinematics type. Use 2 for 2-body, 3 for 3-body, 4 for 4-body!");
		return 1;
	}

	std::cout << std::endl;
	TString msg = Form("Kin%dmc selected!\n",kinX);
	ConsoleColorizer::PrintGreen(msg.Data());

	msg = Form(" Processing physics data file %s\n", argv[2]);
	ConsoleColorizer::PrintGreen(msg.Data());

	msg = Form("Writing to output file %s\n", argv[3]);
	ConsoleColorizer::PrintGreen(msg.Data());

	try{

		SABREsim sim(kinX, argv[2], argv[3]);
		sim.Run();
		return 0;

	} catch(const std::exception& e){
		ConsoleColorizer::PrintRed(Form("Exception: %s\n",e.what()));
		return 1;
	}
	catch(...){
		ConsoleColorizer::PrintRed("Unknown exception occured during simulation!\n");
		return 1;
	}
}
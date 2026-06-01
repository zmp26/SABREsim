#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <vector>
#include <string>

/*
Example command to call this from terminal:
root 'flatten_data.C("../analyze/may31/det/b10ha_7.5MeV_9B_e
x2345keV_p8Be_ex0keV_1mil_tree_mult3_SABREanalyzed8Be.root", "../analyze/may31/det/b10ha_7.5MeV_9B_ex2345keV
_p8Be_ex0keV_1mil_tree_mult3_SABREanalyzed8Be_flattenTest.root")'
*/

void flatten_data(const char* input_filename, const char* output_filename){

	//ROOT::EnableImplicitMT();

	ROOT::RDataFrame df("InvMass_Mult3", input_filename);

	const double simRecEx = 2.345;//this is the "truth" value...i.e. the value passed into kin4mc for the recoil Ex
	std::string simRecExString = std::to_string(simRecEx);
	std::string reconExString = "std::vector<double>{"
								"_012.reconEx - " + simRecExString + ", _021.reconEx - " + simRecExString + ", "
								"_102.reconEx - " + simRecExString + ", _120.reconEx - " + simRecExString + ", "
				  				"_201.reconEx - " + simRecExString + ", _210.reconEx - " + simRecExString + "}";

	auto df_vectors = df
		.Define("reconEx", reconExString)
		.Define("imEx", "std::vector<double>{_012.imEx, _021.imEx, _102.imEx, _120.imEx, _201.imEx, _210.imEx}")
		.Define("delta_ecm1",
			"std::vector<double>{"
				   "_012.ecm1 - _012.exp_ecm1, _021.ecm1 - _021.exp_ecm1,"
				   "_102.ecm1 - _102.exp_ecm1, _120.ecm1 - _120.exp_ecm1,"
				   "_201.ecm1 - _201.exp_ecm1, _210.ecm1 - _210.exp_ecm1"
				   "}"
		)
		.Define("delta_ecm2", 
			"std::vector<double>{"
			"_012.ecm2 - _012.exp_ecm2, _021.ecm2 - _021.exp_ecm2,"
			"_102.ecm2 - _102.exp_ecm2, _120.ecm2 - _120.exp_ecm2,"
			"_201.ecm2 - _201.exp_ecm2, _210.ecm2 - _210.exp_ecm2"
			"}"
		)
		.Define("label", "std::vector<int>{1,0,0,0,0,0}");//for asymmetric decays
		//.Define("label", "std::vector<int>{1,1,0,0,0,0}");//for symmetric decays

	std::vector<std::string> out_cols = {"reconEx", "imEx", "delta_ecm1", "delta_ecm2", "label"};
	df_vectors.Snapshot("SABREtrainTree", output_filename, out_cols);

	std::cout << "Flattening complete with all permutations retained.\nOutput:\t" << output_filename << std::endl;

}
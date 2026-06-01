#include <iostream>
#include <vector>
#include <numeric>
#include "include/onnxruntime_cxx_api.h"

void evaluate_bdt() {

	Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "SABRE_BDT_Inference");
	Ort::SessionOptions session_options;
	session_options.SetIntraOpNumThreads(1);

	const char* model_path = "sabre_bdt.onnx";
	std::cout << "Loading model: " << model_path << std::endl;
	Ort::Session session(env, model_path, session_options);
	//									 {deltaReconEx, imEx, deltaEcm1, deltaEcm2}
	//std::vector<float> input_tensor_values = {2.35f, 12.4f, -0.05f, 0.12f};
	//std::vector<float> input_tensor_values = {-0.022626, 0.0134231, 0.0746820, -6.42e-5};//from 012
	std::vector<float> input_tensor_values = {-0.022626, 0.0134231, 0.0746820, -6.42e-5};//from 021
	std::vector<int64_t> input_node_dims = {1, 4};

	//auto memory_info = Ort::MemoryInfo::CreateCpu(Ort::ArenaAllocator, Ort::MemTypeDefault);
	auto memory_info = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

	Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
		memory_info,
		input_tensor_values.data(),
		input_tensor_values.size(),
		input_node_dims.data(),
		input_node_dims.size()
	);

	const char* input_names[] = {"float_input"};
	const char* output_names[] = {"label", "probabilities"};

	auto output_tensors = session.Run(
		Ort::RunOptions{nullptr},
		input_names,
		&input_tensor, 1,
		output_names, 2
	);


	float* probabilities = output_tensors[1].GetTensorMutableData<float>();

	float prob_true_permutation = probabilities[1];

	std::cout << "--- BDT Evaluation ---" << std::endl;
	std::cout << "Probability of true permutation: " << prob_true_permutation * 100. << "%" << std::endl;
}

int main() {
	evaluate_bdt();
	return 0;
}
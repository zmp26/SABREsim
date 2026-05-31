# SABRE 3-Body Event Permutation BDT Pipeline

This repository contains a complete machine learning pipeline for resolving combinatorics in multiplicity-3 events for the SABRE detector array. The framework flattens 3-body decay permutations, trains an XGBoost Boosted Decision Tree (BDT) classifier to distinguish physical truth from combinatoric background, and exports the model to an optimized C++ executable using the ONNX Runtime API.

---

## Pipeline Architecture and Data Flow

### 1. Data Preparation (`flatten_data.C`)
* **Input:** A ROOT `TTree` (`InvMass_Mult3`) containing multiplicity-3 events analyzed under a specific decay hypothesis. The tree holds structural parameters for all 6 permutations (`_012`, `_021`, `_102`, `_120`, `_201`, `_210`).
* **Operation:** Using **ROOT RDataFrame**, parallel vector columns (`std::vector<double>`) are constructed for the training features. A target vector column is defined as `label = {1, 0, 0, 0, 0, 0}`, explicitly tagging the first permutation (`_012`) as the physical truth and the remaining 5 as combinatoric background.
* **Output:** A flat, multi-threaded snapshot tree (`SABREtrainTree`) containing arrays of fixed dimensions (6 entries per event).

### 2. Model Training and Serialization (`train_bdt.py`)
* **Input:** The flattened ROOT snapshot file.
* **Operation:** Using `uproot`, the vector branches are loaded into memory and completely flattened via `np.concatenate`, removing event boundaries. For example, 181,833 events yield exactly 1,090,998 individual evaluation rows (181,833 events × 6 permutations).
* **Model Fitting:** An **XGBoost Classifier** (`XGBClassifier`) determines the non-linear decision boundaries separating true configurations (Class 1) from background configurations (Class 0) across 4 features: `reconEx`, `imEx`, `delta_ecm1`, and `delta_ecm2`.
* **Serialization:** The trained booster state is exported via `onnxmltools` into a frozen, cross-platform format (`sabre_bdt.onnx`) with a static input footprint (`FloatTensorType`).

### 3. C++ Deployment (`main.cpp` / `bdt_inference`)
* **Input:** The serialized `sabre_bdt.onnx` file and raw runtime feature vectors.
* **Operation:** The **ONNX Runtime (ORT) C++ API** parses the model graph directly into native memory structures. When a 4-feature tensor is passed, ORT traverses the underlying decision trees sequentially, bypassing the Python environment completely.
* **Output:** The network returns a probability distribution over the available classes. `probabilities[1]` provides the exact numerical probability (scaled from `0.0` to `1.0`) that the tested combination matches a true physical decay event.

---

## Environment Replication and Setup

Follow these procedures to configure and compile this framework on any clean Linux or WSL distribution.

### Step 1: Install Python and Training Dependencies
Ensure Miniconda or Anaconda is active, then install the data science and conversion utilities:

```bash
# Update conda to ensure package tracking is modern
conda update -n base -c conda-forge conda -y

# Install Python core libraries and XGBoost
conda install -c conda-forge uproot numpy xgboost -y

# Install ONNX converter libraries via pip
pip install onnxmltools onnxconverter-common

# Execute the training script to generate the serialized model file
python3 train_bdt.py
```

### Step 2: Install Self-Contained C++ ONNX Runtime Headers and Libraries
To eliminate dependency issues related to missing system packages or broken Conda C++ paths, download the official, pre-compiled C/C++ release bundle directly into the local project directory:
```bash
# Fetch the targeted Linux x64 release archive from Microsoft
wget -q --show-progress https://github.com/microsoft/onnxruntime/releases/download/v1.17.1/onnxruntime-linux-x64-1.17.1.tgz

# Extract the 'include' directory containing all necessary .h files
tar -xzvf onnxruntime-linux-x64-1.17.1.tgz --strip-components=1 onnxruntime-linux-x64-1.17.1/include

# Extract the 'lib' directory containing 'libonnxruntime.so'
tar -xzvf onnxruntime-linux-x64-1.17.1.tgz --strip-components=1 onnxruntime-linux-x64-1.17.1/lib

# Remove the raw downloaded archive
rm onnxruntime-linux-x64-1.17.1.tgz
```
The resulting directory structure will conform to the layout below:

├── flatten_data.C
├── main.cpp
├── sabre_bdt.onnx
├── train_bdt.py
├── include/       <-- Local C/C++ Header Files
└── lib/           <-- Local libonnxruntime.so


### Step 3: Compile the C++ Executable
Compile main.cpp by directing g++ to the local include and lib folders. This creates a portable build independent of system environment paths:
```bash
g++ -std=c++17 main.cpp -o bdt_inference \
  -Iinclude \
  -Llib \
  -lonnxruntime \
  -Wl,-rpath,lib
```

### Step 4: Run Inference
Execute the final binary to verify correct runtime execution:
```bash
./bdt_inference
```
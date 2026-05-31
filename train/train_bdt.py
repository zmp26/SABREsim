import uproot
import numpy as np
import xgboost as xgb
import onnxmltools
#from onnxconverter_common.data_types import FloatTensorType
from onnxmltools.convert.common.data_types import FloatTensorType


file = uproot.open("../analyze/may31/det/b10ha_7.5MeV_9B_ex2345keV_p8Be_ex0keV_1mil_tree_mult3_SABREanalyzed8Be_flattenTest.root")
tree = file["SABREtrainTree"]

features_dict = tree.arrays(["reconEx", "imEx", "delta_ecm1", "delta_ecm2"], library="np")
labels_dict = tree.arrays(["label"], library="np")

reconEx_flat = np.concatenate(features_dict["reconEx"])
imEx_flat = np.concatenate(features_dict["imEx"])
delta_ecm1_flat = np.concatenate(features_dict["delta_ecm1"])
delta_ecm2_flat = np.concatenate(features_dict["delta_ecm2"])

X = np.column_stack((reconEx_flat, imEx_flat, delta_ecm1_flat, delta_ecm2_flat)).astype(np.float32)
Y = np.concatenate(labels_dict["label"]).astype(np.int32)

print(f"Training data prepared successfully!")
print(f"Feature matrix X shape: {X.shape} (Total rows = Events * Permutations)")
print(f"Labels array Y  shape: {Y.shape}")

model = xgb.XGBClassifier(
	n_estimators=100,
	max_depth=4,
	learning_rate=0.1,
	objective="binary:logistic"
)
print("Training model...")
model.fit(X,Y)
print("XGBoost training complete!")

initial_type = [('float_input', FloatTensorType([None, X.shape[1]]))]

print("Converting model to ONNX...")
#onnx_model = onnxmltools.convert_xgboost(model.get_booster(), initial_types=initial_type)
onnx_model = onnxmltools.convert_xgboost(model, initial_types=initial_type)

with open("sabre_bdt.onnx", "wb") as f:
    f.write(onnx_model.SerializeToString())

print("Model successfully converted and saved to sabre_bdt.onnx")

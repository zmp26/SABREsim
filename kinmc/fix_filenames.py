import os
import glob

directory = "/mnt/e/SABREsim/det/kin4mc_eff_10B3He4He9B_p_4He"
extension = ".det"

files = glob.glob(os.path.join(directory, f"*{extension}"))

for file_path in files:
	dir_name = os.path.dirname(file_path)
	base_name = os.path.basename(file_path)

	splitted = base_name.split("_")
	print(splitted)
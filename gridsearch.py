import itertools
import subprocess
import tempfile
from pathlib import Path
import re

#CONFIG_PATH = Path("config/SABREsim_laptop.conf")
CONFIG_DIR = Path("config")
LOG_DIR = Path("logs")
SABRESIM_EXE = Path("bin/SABREsim")

CONFIG_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

#helper function to extract gaus value:
def extract_gaus_value(s):
	match = re.search(r'\d+',s)
	if not match:
		return None

	digits = match.group()

	value = float("0."+digits)
	return value

# x_opts = ["fixed", "gaus001", "gaus002", "gaus003", "gaus004", "gaus005"]
# y_opts = ["fixed", "gaus001", "gaus002", "gaus003", "gaus004", "gaus005"]

x_opts = ["fixed", "gaus0005", "gaus001", "gaus0015", "gaus002", "gaus0025"]
y_opts = ["fixed", "gaus0005", "gaus001", "gaus0015"]

for profx in x_opts:

	parx = 0 if profx=="fixed" else extract_gaus_value(profx)

	for profy in y_opts:

		pary = 0 if profy=="fixed" else extract_gaus_value(profy)

		if profx == "fixed" and profy == "fixed":
			profile_str = "fixedpoint"
		elif profx == "fixed":
			profile_str = "fixedxgausy"
		elif profy == "fixed":
			profile_str = "gausxfixedy"
		else:
			profile_str = "gaussian"

		detfilename = f"/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta16282328_{profx}x_{profy}y.det"
		treefilename = f"/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta16282328_{profx}x_{profy}y_tree.root"
		histofilename = f"/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta16282328_{profx}x_{profy}y_histos.root"

		detpath = Path(detfilename)
		treepath = Path(treefilename)
		histopath = Path(histofilename)

		if detpath.exists() and treepath.exists() and histopath.exists():
			print(f"Files detected for profiles x = {profx} and y = {profy} so skipping these...")
			continue

		outputlines = [
					   "#Default SABREsim config file\n\n[General]\ndetmc_version = 2",
					   "infile = /mnt/e/SABREsim/kinmc/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta16282328.out",
					   f"detfile = {detfilename}",
					   f"treefile = {treefilename}",
					   f"histofile = {histofilename}",
					   "\n[TargetLosses]\ntargetLoss_par1 = alpha_in_LiF\ntargetLoss_par2 = 6Li_in_LiF\ntargetLoss_par3 = none\ntargetLoss_par4 = none",
					   "\n[DeadLayerLosses]\ndeadLayerLoss_par1 = alpha_in_Si\ndeadLayerLoss_par2 = 6Li_in_Si\ndeadLayerLoss_par3 = none\ndeadLayerLoss_par4 = none",
					   "\n[Beamspot]",
					   f"profile = {profile_str}",
					   f"parX = {parx}",
					   f"parY = {pary}",
					   "beam_offsetX = 0\nbeam_offsetY = 0",
					   "\n[MetaData]",
					   "reaction = 7Li(3He,4He)6Li\nbeam_energy = 7.5\nrecoil_excitation_energy = 0"
					   ]

		config_file = CONFIG_DIR / f"SABREsim_gridsearch_{profx}_{profy}.conf"
		with tempfile.NamedTemporaryFile(mode="w",suffix=".conf",delete=False) as tmp:
			tmp.write("\n".join(outputlines))
			tmp.flush()
			tmp_path = tmp.name

		log_file = LOG_DIR/f"SABREsim_theta16282328_{profx}x_{profy}y.log"

		with open(log_file, "w") as logf:
			print(f"Running SABREsim with profiles x = {profx}, y = {profy}, and beam parameters x = {parx}, y = {pary}")
			subprocess.run([str(SABRESIM_EXE), tmp_path], stdout=logf, stderr=subprocess.STDOUT)

		Path(tmp_path).unlink()


		# Uncomment below to print out config file names and contents for verification purposes:
		# print(config_file/"\n")
		# for line in outputlines:
		# 	print(line)
		# print("--------------------------------------------------------------")






"""
#Default SABREsim config file

[General]
detmc_version = 2
infile = /mnt/e/SABREsim/kinmc/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta16782278.out
detfile = /mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta16782278_fixedx_fixedy_fixedx_gaus002y.det
treefile = /mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta16782278_fixedx_fixedy_fixedx_gaus002y_tree.root
histofile = /mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta16782278_fixedx_fixedy_fixedx_gaus002y_histos.root

[TargetLosses]
targetLoss_par1 = alpha_in_LiF
targetLoss_par2 = 6Li_in_LiF
targetLoss_par3 = none
targetLoss_par4 = none

[DeadLayerLosses]
deadLayerLoss_par1 = alpha_in_Si
deadLayerLoss_par2 = 6Li_in_Si
deadLayerLoss_par3 = none
deadLayerLoss_par4 = none

[Beamspot]
profile = fixedxgausy
parX = 0.0
parY = 0.002
beam_offsetX = 0
beam_offsetY = 0

[MetaData]
reaction = 7Li(3He,4He)6Li
beam_energy = 7.5
recoil_excitation_energy = 0
"""
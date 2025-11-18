import itertools
import subprocess
import tempfile
from pathlib import Path
import re
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor

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

x_opts = ["gaus001"]
y_opts = ["gaus001"]

#anglestrings = ["188208", "178218", "168228", "158238", "148248"]
#	19.8deg		+/- 1deg  +/- 2deg  +/- 3deg  +/- 4deg  +/- 5deg
anglestrings = [
				#"193203",	# 19.8deg +/- 0.5 deg
				#"188208",	# 19.8deg +/- 1.0 deg
				#"183213",	# 19.8deg +/- 1.5 deg
				"178218",	# 19.8deg +/- 2.0 deg
				#"173223",	# 19.8deg +/- 2.5 deg
				#"168228",	# 19.8deg +/- 3.0 deg
				#"163233",	# 19.8deg +/- 3.5 deg
				#"158238",	# 19.8deg +/- 4.0 deg
				#"153243",	# 19.8deg +/- 4.5 deg
				#"148248"	# 19.8deg +/- 5.0 deg
				]

#kin2mc_7Li3He4He6Ligs_7500keV_theta178218_phi_-0.5_0.5.out
phis = [1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]
#phis = [0.5, 1.0, 1.5, 2.0, 2.5]
phistrings = []
for philow in phis:
	for phiup in phis:
		phistr = f"{-philow}_{phiup}"
		phistrings.append(phistr)
		#print(f"kin2mc_7Li3He4He6Ligs_7500keV_theta178218_phi_{phistr}.out")

#generate the tasks list:
tasks = []
for profx in x_opts:
	for profy in y_opts:
		for angstr in anglestrings:
			for phistring in phistrings:
				tasks.append((profx, profy, angstr, phistring))

def run_simulation(args):
	profx, profy, angstr, phistr = args

	parx = 0 if profx == "fixed" else extract_gaus_value(profx)
	pary = 0 if profy == "fixed" else extract_gaus_value(profy)

	if profx == "fixed" and profy == "fixed":
		profile_str = "fixedpoint"
	elif profx == "fixed":
		profile_str = "fixedxgausy"
	elif profy == "fixed":
		profile_str = "gausxfixedy"
	else:
		profile_str = "gaussian"


	infilename = f"/mnt/e/SABREsim/kinmc/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta{angstr}_phi_{phistr}.out"
	detfilename = f"/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta{angstr}_phi_{phistr}_{profx}x_{profy}y.det"
	treefilename = f"/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta{angstr}_phi_{phistr}_{profx}x_{profy}y_tree.root"
	histofilename = f"/mnt/e/SABREsim/det/kin2mc/kin2mc_7Li3He4He6Ligs_7500keV_theta{angstr}_phi_{phistr}_{profx}x_{profy}y_histos.root"


	outputlines = [
						   "#Default SABREsim config file\n\n[General]\ndetmc_version = 2",
						   f"infile = {infilename}",
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

	with tempfile.NamedTemporaryFile(mode="w", suffix=".conf", delete=False) as tmp:
		tmp.write("\n".join(outputlines))
		tmp.flush()
		tmp_path = tmp.name

	t = datetime.now().strftime("%Y%m%d_%H%M%S")
	log_file = LOG_DIR/f"SABREsim_theta{angstr}_phi_{phistr}.log"

	with open(log_file, "w") as logf:
		print(f"running phistring = {phistr}")
		subprocess.run([str(SABRESIM_EXE), tmp_path], stdout=logf, stderr=subprocess.STDOUT)

	Path(tmp_path).unlink(missing_ok=True)





if __name__ == "__main__":

	start = datetime.now()

	with ProcessPoolExecutor(max_workers=12) as exe:
		exe.map(run_simulation, tasks)

	end = datetime.now()
	elapsed = end - start
	elapsed_sec = elapsed.total_seconds()
	elapsed_min = elapsed_sec/60.

	if elapsed_sec < 60:
		print(f"\nAll SABREsim jobs finished in {elapsed_sec:.2f} seconds")
	else:
		print(f"\nAll SABREsim jobs finished in {elapsed_min:.2f} minutes")
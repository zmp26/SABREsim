import itertools
import subprocess
import tempfile
from pathlib import Path
import re
from datetime import datetime

#CONFIG_PATH = Path("config/SABREsim_laptop.conf")
CONFIG_DIR = Path("config")
LOG_DIR = Path("logs")
SABRESIM_EXE = Path("bin/SABREsim")

CONFIG_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

CONFIG_FILENAME = CONFIG_DIR / "autoSABRE.conf"

def run_simulation(recoilEx_keV):

	recoilEx_MeV = recoilEx_keV/1000.

	datestring = "may16"

	infilename = f"/mnt/e/SABREsim/analyze/{datestring}/kinmc/b10ha_7.5MeV_9B_ex{recoilEx_keV}keV_a5Li_ex0keV_1000k.out"
	detfilename = f"/mnt/e/SABREsim/analyze/{datestring}/det/b10ha_7.5MeV_9B_ex{recoilEx_keV}keV_a5Li_ex0keV_1000k.det"
	treefilename = f"/mnt/e/SABREsim/analyze/{datestring}/det/b10ha_7.5MeV_9B_ex{recoilEx_keV}keV_a5Li_ex0keV_1000k_tree.root"
	histofilename = f"/mnt/e/SABREsim/analyze/{datestring}/det/b10ha_7.5MeV_9B_ex{recoilEx_keV}keV_a5Li_ex0keV_1000k_histos.root"

	outputlines = [
						"[General]\ndetmc_version = 4",
						f"infile = {infilename}",
						f"detfile = {detfilename}",
						f"treefile = {treefilename}",
						f"histofile = {histofilename}",
						f"\n[Kinematics]\nbeam = 3 He\ntarget = 10 B\nrecoil = 9 B\nejectile = 4 He\nbeam_energy = 7.5\nrecoil_excitation_energy = {recoilEx_MeV}",
						"\n[TargetLosses]\ntargetLoss_par1 = none\ntargetLoss_par2 = none\ntargetLoss_par3 = none\ntargetLoss_par4 = none",
						"\n[TargetAngularStraggling]\nenableStraggle1 = false\nenableStraggle2 = false\nenableStraggle3 = false\nenableStraggle4 = false",
						"\n[DeadLayerLosses]\ndeadLayerLoss_par1 = none\ndeadLayerLoss_par2 = none\ndeadLayerLoss_par3 = none\ndeadLayerLoss_par4 = none\n",
						"\n[Beamspot]\nprofile = gaussian\nparX = 0.0005\nparY = 0.0005\nbeam_offsetX = 0\nbeam_offsetY = 0",
						"\n[SPS]\n Coincidence = true\nThetaMin = 18.27\nThetaMax = 21.73\nPhiMin = -2.19\nPhiMax = 2.19\nSigmaE = 0.015\nSigmaTheta = 0.5\nSigmaPhi = 0.5\nApertureDist = 0.175"
				  ]

	with open(CONFIG_FILENAME, "w") as f:
		for line in outputlines:
			f.write(line + "\n")

	print(f"Running SABREsim for {recoilEx_keV} keV...")

	subprocess.run([str(SABRESIM_EXE), str(CONFIG_FILENAME)], check=True)




if __name__ == "__main__":

	start = datetime.now()

	#recoil_values = range(0, 501, 50)
	recoil_values = range(1500, 2501, 50)

	for keV in recoil_values:
		run_simulation(keV)


	end = datetime.now()
	elapsed = end - start
	elapsed_sec = elapsed.total_seconds()
	elapsed_min = elapsed_sec/60.

	if elapsed_sec < 60:
		print(f"\nAll SABREsim jobs finished in {elapsed_sec:.2f} seconds")
	else:
		print(f"\nAll SABREsim jobs finished in {elapsed_min:.2f} minutes")
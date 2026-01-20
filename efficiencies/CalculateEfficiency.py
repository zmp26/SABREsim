#CalculateEfficiency.py
#purpose: 
#	Generates an input file for kin2/3/4mc and pass it as input. A SABREsim config file is generated w/ path of kin2/3/4mc output as the input file.
#   Runs SABREsim w/ the new config file and keeps track of the output (specifically, the number of 1/2/3 particle hits). Does this for a range
#   of excitation energies and generates a txt file w/ the results and a png of the resulting curve.
#
#   Note that all kin2/3/4mc output files are saved, but the input files are overwritten.
#   Note also that all SABREsim output files are saved, but the input config file is overwritten.
#	Logs of the processes can be found in the logs/ subdirectory.
#	
#Result:
#	this code generates a txt file and a plot of the efficiency across the specified energy range for each case
#
#How to run:
#
#	in this directory:
#		python3 CalculateEfficiency.py X
#
#	where X = (2,3,4)
#
#
#
# The paths should be updated in the code before running. This feels easier than passing them as command line arguments...but maybe not

import sys
import os
import subprocess
import re
import matplotlib.pyplot as plt
import time
import logging


'''
Path Constants (per run/resonance):
'''

#kinmc
KINMC_INPUT_FILEPATH = "/mnt/e/SABREsim/efficiencies/kinmc.in"
KINMC_OUTPUT_DIR = "/mnt/e/SABREsim/efficiencies/kinmc_output/"

#sabresim config file:
SABRESIM_CONFIG_FILEPATH = "/mnt/e/SABREsim/config/EffCalc.conf"
SABRESIM_CONFIG_DETDIR = "/mnt/e/SABREsim/efficiencies/det/"
SABRESIM_CONFIG_TREEDIR = "/mnt/e/SABREsim/efficiencies/tree/"
SABRESIM_CONFIG_HISTODIR = "/mnt/e/SABREsim/efficiencies/histos/"
#sabresim logging:
SABRESIM_LOGDIR = "/mnt/e/SABREsim/logs/"
#sabresim executable:
SABRESIM_EXE = "/mnt/e/SABREsim/bin/SABREsim"

#efficiency calculation log - automation log
LOG_FILEPATH = "/mnt/e/SABREsim/efficiencies/logs/"

EFFICIENCIES_FILEPATH = "/mnt/e/SABREsim/efficiencies/efficiencies.txt"
EFFICIENCIES_PLOT_FILEPATH = "/mnt/e/SABREsim/efficiencies/efficiencies.png"

'''
Resonance/reaction specific SABREsim config data:
'''
#targetlosses:
TARGETLOSS_PAR1 = "alpha_in_LiF"
TARGETLOSS_PAR2 = "6Li_in_LiF"
TARGETLOSS_PAR3 = "alpha_in_LiF"
TARGETLOSS_PAR4 = "deuteron_in_LiF"

#target angular straggling:
ENABLESTRAGGLE1 = False
ENABLESTRAGGLE2 = False
ENABLESTRAGGLE3 = False
ENABLESTRAGGLE4 = False

STRAGGLE1 = "none"
STRAGGLE2 = "none"
STRAGGLE3 = "none"
STRAGGLE4 = "none"

#dead layer losses
DEADLAYERLOSS_PAR1 = "alpha_in_Si"
DEADLAYERLOSS_PAR2 = "6Li_in_Si"
DEADLAYERLOSS_PAR3 = "alpha_in_Si"
DEADLAYERLOSS_PAR4 = "deuteron_in_Si"

#beam spot
PROFILE = "gaussian"
PARX = "0.0005"
PARY = "0.0005"
BEAMOFFSETX = "0"
BEAMOFFSETY = "0"

#meta data - we actually update this in kin2/3/4mc

'''
Binning:

'''

ENERGY_START_KEV = 1500
ENERGY_STOP_KEV = 1510
ENERGY_BIN_KEV = 5 #bin size





#create the globally accessible log file here:
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def extractEfficienciesFromSTDOUT_kin2mc(stdout):
	eff = {'ej':None, 'rec':None, 'both':None}
	lines = stdout.splitlines()
	for line in lines:
		low = line.lower()
		if "only ejectile" in low:
			match = re.search(r"\(([\d.]+)%", low)
			if match:
				eff['ej'] = float(match.group(1))
		elif "only recoil" in low:
			match = re.search(r"\(([\d.]+)%", low)
			if match:
				eff['rec'] = float(match.group(1))
		elif "both" in low:
			match = re.search(r"\(([\d.]+)%", low)
			if match:
				eff['both'] = float(match.group(1))
	return eff

def kin2mc():
	logger.info("kin2mc() started")

	FIXED_OUTPUT_NAME = "kin2mc.out" #this is the default name of kin2mc output which we move and rename - don't change this unless kin2mc changes!

	#prepare information needed for kin2mc input file (all keystrokes and inputs otherwise manually entered when using kin2mc)
	energy_start_keV = float(ENERGY_START_KEV)
	energy_stop_keV = float(ENERGY_STOP_KEV)
	energy_step_keV = ENERGY_BIN_KEV

	target = "7Li" #UPDATE THIS PER CURVE
	beam = "3He" #UPDATE THIS PER CURVE
	ejectile = "4He" #UPDATE THIS PER CURVE
	recoil = "6Li" #UPDATE THIS PER CURVE

	beamenergy = 7.5 #UPDATE THIS PER CURVE

	reaction = target + "(" + beam + "," + ejectile + ")" + recoil #this is for kin2mc input purposes
	reaction_nospecchars = target + beam + ejectile + recoil #this is for filename purposes

	exc1_MeV = 0. #ejectile excitation energy - update this per curve if necessary
	exc2_keV = energy_start_keV #recoil excitation energy - update this per cycle in main loop
	thetamin_deg = 0.
	thetamax_deg = 180.
	phimin_deg = 0.
	phimax_deg = 360.
	use_lookup = 0
	num_events = 1000000

	#main loop here:
	results = []
	counter = 0
	while exc2_keV <= energy_stop_keV:
		logger.info(f"readying to run kin2mc for ExE = {exc2_keV} keV")
		input_lines = [
			"1",
			reaction,
			"3",
			str(beamenergy),
			"4",
			str(exc1_MeV),
			str(exc2_keV/1000.),
			"5",
			str(thetamin_deg),
			str(thetamax_deg),
			str(phimin_deg),
			str(phimax_deg),
			str(use_lookup),
			"6",
			str(num_events),
			"9"
		]

		in_filename = KINMC_INPUT_FILEPATH
		out_filename = os.path.join(KINMC_OUTPUT_DIR,f"kin2mc_{reaction_nospecchars}_ExE{int(round(exc2_keV))}keV.out")

		with open(in_filename, "w") as f:
			for line in input_lines:
				f.write(line + "\n")

		logger.info(f"passing kin2mc input file to kin2mc for ExE = {exc2_keV} keV")
		#subprocess.run(["/mnt/e/kinematics/kin2mc/kin2mc_legacy/src/kin2mc"], stdin=open(in_filename),check=True)
		with open(in_filename) as fin:
			subprocess.run(["/mnt/e/kinematics/kin2mc/kin2mc_legacy/src/kin2mc"], stdin=fin, check=True)
		os.replace(FIXED_OUTPUT_NAME, out_filename)
		logger.info(f"finished running kin2mc for ExE = {exc2_keV} keV")


		logger.info(f"readying SABREsim config file for ExE = {exc2_keV} keV")
		infilename = out_filename #sabresim in file is the outfile of kin2mc from above!
		detfilename = os.path.join(SABRESIM_CONFIG_DETDIR, f"kin2mc_{reaction_nospecchars}_ExE{round(exc2_keV)}keV.det")
		treefilename = os.path.join(SABRESIM_CONFIG_TREEDIR, f"kin2mc_{reaction_nospecchars}_ExE{round(exc2_keV)}keV_tree.root")
		histofilename = os.path.join(SABRESIM_CONFIG_HISTODIR, f"kin2mc_{reaction_nospecchars}_ExE{round(exc2_keV)}keV_histos.root")

		es1 = "true" if ENABLESTRAGGLE1 else "false"
		es2 = "true" if ENABLESTRAGGLE2 else "false"
		es3 = "true" if ENABLESTRAGGLE3 else "false"
		es4 = "true" if ENABLESTRAGGLE4 else "false"

		outputlines = [
			"#Default SABREsim config file\n\n[General]\ndetmc_version = 2",
			f"\ninfile = {infilename}",
			f"\ndetfile = {detfilename}",
			f"\ntreefile = {treefilename}",
			f"\nhistofile = {histofilename}",
			"\n\n[TargetLosses]",
			f"\ntargetLoss_par1 = {TARGETLOSS_PAR1}",
			f"\ntargetLoss_par2 = {TARGETLOSS_PAR2}",
			f"\ntargetLoss_par3 = {TARGETLOSS_PAR3}",
			f"\ntargetLoss_par4 = {TARGETLOSS_PAR4}",
			"\n\n[TargetAngularStraggling]",
			f"\nenableStraggle1 = {es1}",
			f"\nstraggle1 = {STRAGGLE1}",
			f"\nenableStraggle2 = {es2}",
			f"\nstraggle2 = {STRAGGLE2}",
			f"\nenableStraggle3 = {es3}",
			f"\nstraggle3 = {STRAGGLE3}",
			f"\nenableStraggle4 = {es4}",
			f"\nstraggle4 = {STRAGGLE4}",
			"\n\n[DeadLayerLosses]",
			f"\ndeadLayerLoss_par1 = {DEADLAYERLOSS_PAR1}",
			f"\ndeadLayerLoss_par2 = {DEADLAYERLOSS_PAR2}",
			f"\ndeadLayerLoss_par3 = {DEADLAYERLOSS_PAR3}",
			f"\ndeadLayerLoss_par4 = {DEADLAYERLOSS_PAR4}",
			"\n\n[Beamspot]",
			f"\nprofile = {PROFILE}",
			f"\nparX = {PARX}",
			f"\nparY = {PARY}",
			f"\nbeam_offsetX = {BEAMOFFSETX}",
			f"\nbeam_offsetY = {BEAMOFFSETY}",
			"\n\n[MetaData]",
			f"\nreaction = {reaction}",
			f"\nbeam_energy = {beamenergy}",
			f"\nrecoil_excitation_energy = {exc2_keV/1000.}",
			"\n"
		]

		with open(SABRESIM_CONFIG_FILEPATH, "w") as configfile:
			for line in outputlines:
				configfile.write(line)
		logger.info(f"wrote sabre config file for ExE = {exc2_keV} keV, preparing to run")

		sabresimlogfile = os.path.join(SABRESIM_LOGDIR,f"templog{int(exc2_keV)}keV.txt")
		command = [SABRESIM_EXE, SABRESIM_CONFIG_FILEPATH]
		logger.info(f"Running SABRE command: {command}")
		result = subprocess.run(command, capture_output=True, text=True)

		if result.returncode != 0:
			logger.warning(f"Error running SABREsim for ExE = {exc2_keV} keV")
			exc2_keV += energy_step_keV
			counter += 1
			continue

		with open(sabresimlogfile,"w") as logf:
			if result.stdout:
				logf.write("=== STDOUT ===\n" + result.stdout)
			if result.stderr:
				logf.write("=== STDERR ===\n" + result.stderr)
		logger.info(f"wrote SABRE stdout for ExE = {exc2_keV/1000.} to {sabresimlogfile}")


		efficiencies = extractEfficienciesFromSTDOUT_kin2mc(result.stdout)
		results.append(
			{
				'energy':exc2_keV/1000.,
				'ej':efficiencies['ej'],
				'rec':efficiencies['rec'],
				'both':efficiencies['both']
			}
		)
		logger.info(f"efficiencies extracted and stored for ExE = {exc2_keV} keV")


		exc2_keV += energy_step_keV
		logger.info(f"incremented ExE to {exc2_keV} keV")
		counter += 1

	logger.info(f"finished loop in {counter} iterations")
	#results.sort(key=lambda x: x['energy'])

	energies = [r['energy']*1000. for r in results]
	ej_vals = [r['ej'] for r in results]
	rec_vals = [r['rec'] for r in results]
	both_vals = [r['both'] for r in results]

	plt.plot(energies, ej_vals, label="Ejectile", marker="o")
	plt.plot(energies, rec_vals, label="Recoil", marker="o")
	plt.plot(energies, both_vals, label="Both", marker="o")

	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Efficiency (%)")
	plt.title("Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig(EFFICIENCIES_PLOT_FILEPATH)

	logger.info(f"finished plotting, saved to {EFFICIENCIES_PLOT_FILEPATH}")

	with open(EFFICIENCIES_FILEPATH, "w") as f:
		f.write("ExE(keV)\tej\trec\tboth\n")
		for x in range(len(energies)):
			f.write(f"{energies[x]}\t{ej_vals[x]}\t{rec_vals[x]}\t{both_vals[x]}\n")

	logger.info(f"finished saving efficiencies to {EFFICIENCIES_FILEPATH}")
	logger.info(f"Successfully finished calculating kin2mc efficiencies for {reaction}")


def kin3mc():
	logger.info("kin3mc() started")

	FIXED_OUTPUT_NAME = "kin3mc_bw.out" #this is the default name of kin2mc output which we move and rename - don't change this unless kin3mc changes!

	#prepare information needed for kin3mc input file (all keystrokes and inputs otherwise manually entered when using kin2mc)
	energy_start_keV = float(ENERGY_START_KEV)
	energy_stop_keV = float(ENERGY_STOP_KEV)
	energy_step_keV = ENERGY_BIN_KEV

	target = "7Li"#UPDATE THIS PER CURVE
	beam = "3He"#UPDATE THIS PER CURVE
	ejectile = "4He"#UPDATE THIS PER CURVE
	recoil = "6Li"#UPDATE THIS PER CURVE

	reaction = target + "(" + beam + "," + ejectile + ")" + recoil #this is for kin2mc input purposes
	reaction_nospecchars = target + beam + ejectile + recoil #this is for filename purposes

	decay_particle = "4He"#UPDATE THIS PER CURVE
	decay_half_width_MeV = 0.012#UPDATE THIS PER CURVE


	beamenergy = 7.5#UPDATE THIS PER CURVE
	exc1_MeV = 0.
	exc2_keV = ENERGY_START_KEV
	exc3_MeV = 0.
	exc4_MeV = 0.
	thetamin_deg = 0.
	thetamax_deg = 180.
	phimin_deg = 0.
	phimax_deg = 360.
	use_lookup = 0
	num_events = 100000

	anything = "z" #this is just to continue after kin3mc finishes as it requires a carriage return to finish completely

	results = []
	counter = 0
	while exc2_keV <= energy_stop_keV:
		logger.info(f"readying to run kin3mc for ExE = {exc2_keV} keV")
		input_lines = [
			"1",
			reaction,
			decay_particle,
			"3",
			str(beamenergy),	
			"4",
			str(exc1_MeV),
			str(exc2_keV/1000.),
			str(exc3_MeV),
			str(exc4_MeV),
			str(decay_half_width_MeV),
			"5",
			str(thetamin_deg),
			str(thetamax_deg),
			str(phimin_deg),
			str(phimax_deg),
			"6",
			str(use_lookup),
			str(num_events),
			anything,
			"9"
		]

		in_filename = KINMC_INPUT_FILEPATH
		out_filename = os.path.join(KINMC_OUTPUT_DIR, f"kin3mc_{reaction_nospecchars}_ExE{int(round(exc2_keV))}keV.out")

		with open(in_filename, "w") as f:
			for line in input_lines:
				f.write(line + "\n")


		logger.info(f"passing kin3mc input file to kin3mc for ExE = {exc2_keV} keV")
		with open(in_filename) as fin:
			subprocess.run(["/mnt/e/kinematics/kin3mc/kin3mc_legacy/src/kin3mc"], stdin=fin, check=True)
		os.replace(FIXED_OUTPUT_NAME, out_filename)
		logger.info(f"finished running kin3mc for ExE = {exc2_keV} keV")

		logger.info(f"readying SABREsim config file for ExE = {exc2_keV} keV")
		infilename = out_filename
		detfilename = os.path.join(SABRESIM_CONFIG_DETDIR, f"kin3mc_{reaction_nospecchars}_ExE{round(exc2_keV)}keV.det")
		treefilename = os.path.join(SABRESIM_CONFIG_TREEDIR, f"kin3mc_{reaction_nospecchars}_ExE{round(exc2_keV)}keV_tree.root")
		histofilename = os.path.join(SABRESIM_CONFIG_HISTODIR, f"kin3mc_{reaction_nospecchars}_ExE{round(exc2_keV)}keV_histos.root")

		es1 = "true" if ENABLESTRAGGLE1 else "false"
		es2 = "true" if ENABLESTRAGGLE2 else "false"
		es3 = "true" if ENABLESTRAGGLE3 else "false"
		es4 = "true" if ENABLESTRAGGLE4 else "false"

		outputlines = [
			"#Default SABREsim config file\n\n[General]\ndetmc_version = 3",
			f"\ninfile = {infilename}",
			f"\ndetfile = {detfilename}",
			f"\ntreefile = {treefilename}",
			f"\nhistofile = {histofilename}",
			"\n\n[TargetLosses]",
			f"\ntargetLoss_par1 = {TARGETLOSS_PAR1}",
			f"\ntargetLoss_par2 = {TARGETLOSS_PAR2}",
			f"\ntargetLoss_par3 = {TARGETLOSS_PAR3}",
			f"\ntargetLoss_par4 = {TARGETLOSS_PAR4}",
			"\n\n[TargetAngularStraggling]",
			f"\nenableStraggle1 = {es1}",
			f"\nstraggle1 = {STRAGGLE1}",
			f"\nenableStraggle2 = {es2}",
			f"\nstraggle2 = {STRAGGLE2}",
			f"\nenableStraggle3 = {es3}",
			f"\nstraggle3 = {STRAGGLE3}",
			f"\nenableStraggle4 = {es4}",
			f"\nstraggle4 = {STRAGGLE4}",
			"\n\n[DeadLayerLosses]",
			f"\ndeadLayerLoss_par1 = {DEADLAYERLOSS_PAR1}",
			f"\ndeadLayerLoss_par2 = {DEADLAYERLOSS_PAR2}",
			f"\ndeadLayerLoss_par3 = {DEADLAYERLOSS_PAR3}",
			f"\ndeadLayerLoss_par4 = {DEADLAYERLOSS_PAR4}",
			"\n\n[Beamspot]",
			f"\nprofile = {PROFILE}",
			f"\nparX = {PARX}",
			f"\nparY = {PARY}",
			f"\nbeam_offsetX = {BEAMOFFSETX}",
			f"\nbeam_offsetY = {BEAMOFFSETY}",
			"\n\n[MetaData]",
			f"\nreaction = {reaction}",
			f"\nbeam_energy = {beamenergy}",
			f"\nrecoil_excitation_energy = {exc2_keV/1000.}",
			"\n"
		]

		with open(SABRESIM_CONFIG_FILEPATH, "w") as configfile:
			for line in outputlines:
				configfile.write(line)
		logger.info(f"wrote SABRE config file for ExE = {exc2_keV} keV, preparing to run")

		sabresimlogfile = os.path.join(SABRESIM_LOGDIR,f"templog{int(exc2_keV)}keV.txt")
		command = [SABRESIM_EXE, SABRESIM_CONFIG_FILEPATH]
		logger.info(f"running SABREsim command: {command}")
		result = subprocess.run(command, capture_output=True, text=True)

		if result.returncode != 0:
			logger.warning(f"Error running SABREsim for ExE = {exc2_keV} keV")
			exc2_keV += energy_step_keV
			counter += 1
			continue

		with open(sabresimlogfile,"w") as logf:
			if result.stdout:
				logf.write("=== STDOUT ===\n" + result.stdout)
			if result.stderr:
				logf.write("=== STDERR ===\n" + result.stderr)
		logger.info(f"wrote SABRE stdout for ExE = {exc2_keV/1000.} to {sabresimlogfile}")

		efficiencies = extractEfficienciesFromSTDOUT_kin3mc(result.stdout)
		results.append(
			{
				'energy':exc2_keV/1000.,
				'bu1':efficiencies['bu1'],
				'bu2':efficiencies['bu2'],
				'both':efficiencies['both']
			}
		)
		logger.info(f"efficiencies extracted and stored for ExE = {exc2_keV} keV")


		exc2_keV += energy_step_keV
		logger.info(f"incremented ExE to {exc2_keV} keV")
		counter += 1

	logger.info(f"finished loop in {counter} iterations")
	#results.sort(key=lambda x: x['energy'])

	energies = [r['energy']*1000. for r in results]
	bu1_vals = [r['bu1'] for r in results]
	bu2_vals = [r['bu2'] for r in results]
	both_vals = [r['both'] for r in results]

	plt.plot(energies, ej_vals, label="Ejectile", marker="o")
	plt.plot(energies, bu1_vals, label="Recoil", marker="o")
	plt.plot(energies, bu2_vals, label="Both", marker="o")

	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Efficiency (%)")
	plt.title("Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig(EFFICIENCIES_PLOT_FILEPATH)

	logger.info(f"finished plotting, saved to {EFFICIENCIES_PLOT_FILEPATH}")

	with open(EFFICIENCIES_FILEPATH, "w") as f:
		f.write("ExE(keV)\tbu1\tbu2\tboth\n")
		for x in range(len(energies)):
			f.write(f"{energies[x]}\t{bu1_vals[x]}\t{bu2_vals[x]}\t{both_vals[x]}\n")

	logger.info(f"finished saving efficiencies to {EFFICIENCIES_FILEPATH}")
	logger.info(f"Successfully finished calculating kin3mc efficiencies for {reaction}")

	

def kin4mc():
	logger.info("kin4mc() started")
	

def main():
	numarg = len(sys.argv)
	if(numarg != 2):
		print("Invalid number of arguments! See CalculateEfficiency.py for more info!")
		return

	#file:
	fh = logging.FileHandler('logs/CalculateEfficiency_log.txt')
	fh.setLevel(logging.INFO)

	#console:
	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)

	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)

	logger.addHandler(fh)
	logger.addHandler(ch)

	#use logging via:
	#	logger.info("logging info to file and console")
	#	logger.warning("")
	#	logger.debug("")


	if(sys.argv[1] == '2'):
		print("kin2mc option chosen")
		kin2mc()

	elif(sys.argv[1] == '3'):
		print("kin3mc option chosen")
		kin3mc()

	elif(sys.argv[2] == '4'):
		print("kin4mc option chosen")
		kin4mc()

	else:
		print("Invalid kinXmc argument. Expected 2, 3, or 4 but got " + str(sys.argv[1]))
		return


if __name__ == "__main__":
	main()
#calculateEfficiency.py
#purpose: 
#	to take all of the kinXmc output files generated in the generateInputFilesAndRun.py code
#	and pass them to SABREsim and extract the calculated efficiency for all three cases (detect only bu1, bu2, or both)
#	
#Result:
#	this code generates a txt file and a plot of the efficiency across the specified energy range for each case
#
#How to run:
#
#	in this directory:
#		python3 calculateEfficiency.py X path/to/kinXmc/files/ path/to/SABREsim/outputs/
#
#	where X = (2,3,4)
#	and the paths are paths to directories containing the files, not a files itself!
#

import sys
import os
import subprocess
import re
import matplotlib.pyplot as plt
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

def extractEnergyFromFilename(filename):
	match = re.search(r'at(\d+)keV',filename)
	if match:
		return int(match.group(1))
	return None

def extractEnergyFromFilename_kin2mc(filename):
	# match = re.search(r'kin3mc_(\d+)keV',filename)
	# if match:
	# 	energy_keV = int(match.group(1))
	# 	return energy_keV
	# return None
	matches = re.findall(r'(\d+)keV',filename)
	if len(matches) >= 2:
		resonance_energy_keV = int(matches[1])
		return resonance_energy_keV
	return None

def extractEnergyFromFilename_kin4mc(filename):
	matches = re.findall(r'(\d+)keV',filename)
	if len(matches) >= 2:
		resonance_energy_keV = int(matches[1])
		return resonance_energy_keV
	return None

def extractEfficienciesFromSTDOUT_kin2mc(stdout):
	eff = {'ej':None, 'rec':None, 'both':None}
	lines = stdout.splitlines()
	for line in lines:
		if "only ejectile" in line:
			match = re.search(r"\(([\d.]+)%", line)
			if match:
				eff['ej'] = float(match.group(1))
		elif "only recoil" in line:
			match = re.search(r"\(([\d.]+)%", line)
			if match:
				eff['rec'] = float(match.group(1))
		elif "both" in line:
			match = re.search(r"\(([\d.]+)%", line)
			if match:
				eff['both'] = float(match.group(1))
	return eff

def extractEfficienciesFromSTDOUT_kin3mc(stdout):
	eff = {'bu1':None,'bu2':None,'both':None}
	lines = stdout.splitlines()
	for line in lines:
		if "Only bu1" in line:
			match = re.search(r"\(([\d.]+)%",line)
			if match:
				eff['bu1'] = float(match.group(1))
		elif "Only bu2" in line:
			match = re.search(r"\(([\d.]+)%",line)
			if match:
				eff['bu2'] = float(match.group(1))
		elif "Both bu1 & bu2" in line:
			match = re.search(r"\(([\d.]+)%",line)
			if match:
				eff['both'] = float(match.group(1))
	return eff

def extractEfficienciesFromSTDOUT_kin4mc(stdout):
	eff = {
		   'bu1':None,
		   'bu2':None,
		   'bu3':None,
		   'bu1_bu2':None,
		   'bu2_bu3':None,
		   'bu1_bu3':None,
		   'all':None
	}

	total_events_match = re.search(r"Processed (\d+) kin4mc events",stdout);
	if not total_events_match:
		return eff

	total_events = int(total_events_match.group(1))
	lines = stdout.splitlines()

	for line in lines:
		if "Only bu1:" in line:
			match = re.search(r":\s*(\d+)",line)
			if match:
				count = int(match.group(1))
				eff['bu1'] = (count / total_events) * 100.

		elif "Only bu2:" in line:
			match = re.search(r":\s*(\d+)",line)
			if match:
				count = int(match.group(1))
				eff['bu2'] = (count / total_events) * 100.

		elif "Only bu3:" in line:
			match = re.search(r":\s*(\d+)",line)
			if match:
				count = int(match.group(1))
				eff['bu3'] = (count / total_events) * 100.

		elif "Only bu1 & bu2:" in line:
			match = re.search(r":\s*(\d+)",line)
			if match:
				count = int(match.group(1))
				eff['bu1_bu2'] = (count / total_events) * 100.
		elif "Only bu2 & bu3:" in line:
			match = re.search(r":\s*(\d+)",line)
			if match:
				count = int(match.group(1))
				eff['bu2_bu3'] = (count / total_events) * 100.

		elif "Only bu1 & bu3:" in line:
			match = re.search(r":\s*(\d+)",line)
			if match:
				count = int(match.group(1))
				eff['bu1_bu3'] = (count / total_events) * 100.

		elif "Only bu1, bu2 & bu3:" in line:
			match = re.search(r":\s*(\d+)",line)
			if match:
				count = int(match.group(1))
				eff['all'] = (count / total_events) * 100.


	return eff

def run_SABREsim_kin4mc(fn, input_dir, output_dir):
	"""Run SABREsim for single file and extract results"""
	input_file = os.path.join(input_dir,fn)
	output_file = os.path.join(output_dir,fn.replace(".out",".det"))

	print(f"\n----------------------------------------------------------------------------------\n")
	print(f"input file = {input_file}")
	print(f"output_file = {output_file}")

	command = ["bin/SABREsim", "4", input_file, output_file]
	print(f"Running command: {' '.join(command)}")

	result = subprocess.run(command, capture_output=True, text=True)
	if(result.returncode != 0):
		print(f"Error running SABREsim fore {fn}")
		print("STDERR:". result.stderr)
		return None

	energy = extractEnergyFromFilename_kin4mc(fn)
	efficiencies = extractEfficienciesFromSTDOUT_kin4mc(result.stdout)

	if energy is None:
		print(f"returning None for {fn}")
		return None

	print(
			f"bu1 eff = {efficiencies['bu1']:.4f}\n"
			f"bu2 eff = {efficiencies['bu2']:.4f}\n"
			f"bu3 eff = {efficiencies['bu3']:.4f}\n"
			f"bu1_bu2 eff = {efficiencies['bu1_bu2']:.4f}\n"
			f"bu2_bu3 eff = {efficiencies['bu2_bu3']:.4f}\n"
			f"bu1_bu3 eff = {efficiencies['bu1_bu3']:.4f}\n"
			f"all eff = {efficiencies['all']:.4f}\n"
		)

	return {
		'energy':energy,
		'bu1':efficiencies['bu1'],
		'bu2':efficiencies['bu2'],
		'bu3':efficiencies['bu3'],
		'bu1_bu2':efficiencies['bu1_bu2'],
		'bu2_bu3':efficiencies['bu2_bu3'],
		'bu1_bu3':efficiencies['bu1_bu3'],
		'all':efficiencies['all'],
	}





def kin2mc():
	print("kin2mc() running...")
	start_time = time.time()

	input_dir = os.path.abspath(sys.argv[2])
	output_dir = os.path.abspath(sys.argv[3])

	os.makedirs(output_dir,exist_ok=True)

	filenames = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir,f))]

	results = []

	count = 0

	for fn in filenames:
		input_file = os.path.join(input_dir,fn)
		output_file = os.path.join(output_dir,fn.replace(".out",".det"))
		#print("test1")
		print("\n----------------------------------------------------------------------------------\n")

		print(f"input_file = {input_file}")
		print(f"output_file = {output_file}")

		command = [
		"../src/SABREsim",
		"2",
		input_file,
		output_file
		]

		print(f"Running command: {' '.join(command)}")

		result = subprocess.run(command, capture_output=True, text=True)

		if result.returncode != 0:
			print(f"Error running SABREsim for {fn}")
			print("STDERR: ", result.stderr)
			continue

		energy = extractEnergyFromFilename(fn)
		efficiencies = extractEfficienciesFromSTDOUT_kin2mc(result.stdout)

		if energy is not None:
			results.append({
				'energy':energy,
				'ej':efficiencies['ej'],
				'rec':efficiencies['rec'],
				'both':efficiencies['both']
				})
			count += 1
		#print("test")
	results.sort(key=lambda x: x['energy'])

	#plot!
	energies = [r['energy'] for r in results]
	ej_vals = [r['ej'] for r in results]
	rec_vals = [r['rec'] for r in results]
	both_vals = [r['both'] for r in results]

	plt.plot(energies, ej_vals, label="Ejectile",marker="o")
	plt.plot(energies, rec_vals, label="Recoil",marker="o")
	plt.plot(energies, both_vals, label="Both",marker="o")

	plt.xlabel("Beam Energy (keV)")
	plt.ylabel("Efficiency (%)")
	plt.title("Detection Efficiency vs Beam Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig(os.path.join(output_dir,"efficiency_plot.png"))

	stop_time = time.time()
	elapsed_time = stop_time - start_time
	print(f"\n\nFinished running {count} files through SABREsim in {elapsed_time:.2f} seconds")

	plt.show()

def kin3mc():
	print("kin3mc() running...")
	#input_dir = sys.argv[2]
	start_time = time.time()

	input_dir = os.path.abspath(sys.argv[2])
	output_dir = os.path.abspath(sys.argv[3])

	os.makedirs(output_dir,exist_ok=True)

	filenames = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

	results = []

	count = 0

	fn_counter = 0
	filenames_size = len(filenames)
	for fn in filenames:
		input_file = os.path.join(input_dir, fn)
		output_file = os.path.join(output_dir, fn.replace(".out",".det"))

		print("\n----------------------------------------------------------------------------------\n")

		print(f"input_file = {input_file}")
		print(f"output_file = {output_file}")

		command = [
		"../src/SABREsim",
		"3",
		input_file,
		output_file
		]

		print(f"Running command: {' '.join(command)}")

		result = subprocess.run(command, capture_output=True, text=True)

		if result.returncode != 0:
			print(f"Error running SABREsim for {fn}")
			print("STDERR:", result.stderr)
			continue
		
		energy = extractEnergyFromFilename_kin2mc(fn)
		efficiencies = extractEfficienciesFromSTDOUT_kin3mc(result.stdout)

		if energy is not None:
			results.append({
				'energy':energy,
				'bu1':efficiencies['bu1'],
				'bu2':efficiencies['bu2'],
				'both':efficiencies['both']
				})
			count += 1

		#breaka
		print(f"bu1 eff = {efficiencies['bu1']}\nbu2 eff = {efficiencies['bu2']}\nboth eff = {efficiencies['both']}")
		fn_counter += 1
		print(f"Finished {fn_counter} of {filenames_size}\t({fn_counter*100./filenames_size}%)")
		#break

	#sort by energy
	# results.sort(key=lambda x: x['energy'])

	# #plot!
	# energies = [r['energy'] for r in results]
	# bu1_vals = [r['bu1'] for r in results]
	# bu2_vals = [r['bu2'] for r in results]
	# both_vals = [r['both'] for r in results]
	# bu1_both_vals = [(r['bu1'] + r['both']) for r in results]
	# bu2_both_vals = [(r['bu2'] + r['both']) for r in results]

	# #normalized:
	# maxval = max(bu1_vals)
	# bu1_vals_norm = [val/maxval for val in bu1_vals]

	# maxval = max(bu2_vals)
	# bu2_vals_norm = [val/maxval for val in bu2_vals]

	# maxval = max(both_vals)
	# both_vals_norm = [val/maxval for val in both_vals]

	# #plt.plot(energies, bu1_vals, label="Only bu1", marker='o')
	# #plt.plot(energies, bu2_vals, label="Only bu2", marker='o')
	# #plt.plot(energies, both_vals, label="Both", marker='o')
	# #plt.plot(energies, bu1_both_vals, label="Only bu1 + both", marker="o")
	# #plt.plot(energies, bu2_both_vals, label="Only bu2 + both", marker="o")
	# plt.plot(energies, bu1_vals_norm, label="Only bu1", marker="o")
	# plt.plot(energies, bu2_vals_norm, label="Only bu2", marker="o")
	# plt.plot(energies, both_vals_norm, label="Both bu1 and bu2, normalized", marker="o")
	# plt.xlabel("Beam Energy (keV)")
	# plt.ylabel("Relative Efficiency (%)")
	# plt.title("Detection Efficiency vs Beam Energy")
	# plt.legend()
	# plt.grid(True)
	# plt.tight_layout()
	# plt.savefig(os.path.join(output_dir,"efficiency_plot.png"))

	# txt_output_path = os.path.join(output_dir, "efficiency_plot_data.txt")
	# with open(txt_output_path, "w") as f:
	# 	f.write("Energy(keV)\tbu1_eff\tbu1_eff_norm\tbu2_eff\tbu2_eff_norm\tboth_eff\tboth_eff_norm\n")
	# 	for i in range(len(energies)):
	# 		f.write(f"{energies[i]}\t{bu1_vals[i]:.4f}\t{bu1_vals_norm[i]:.4f}\t{bu2_vals[i]:.4f}\t{bu2_vals_norm[i]:.4f}\t{both_vals[i]:.4f}\t{both_vals_norm[i]:.4f}\n")
	# print(f"\nSaved efficiency plot data to: {txt_output_path}\n")

	# stop_time = time.time()
	# elapsed_time = stop_time - start_time
	# print(f"\n\nFinished running {count} files through SABREsim in {elapsed_time:.2f} seconds")

	# plt.show()



# def kin4mc():
# 	print("kin4mc() running...")

# 	start_time = time.time()

# 	input_dir = os.path.abspath(sys.argv[2])
# 	output_dir = os.path.abspath(sys.argv[3])

# 	os.makedirs(output_dir,exist_ok=True)

# 	filenames = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir,f))]

# 	results = []

# 	count = 0

# 	fn_counter = 0
# 	filenames_size = len(filenames)
# 	for fn in filenames:
# 		input_file = os.path.join(input_dir,fn)
# 		output_file = os.path.join(output_dir,fn.replace(".out",".det"))

# 		print("\n----------------------------------------------------------------------------------\n")

# 		print(f"input_file = {input_file}")
# 		print(f"output_file = {output_file}")

# 		command = [
# 		"bin/SABREsim",
# 		"4",
# 		input_file,
# 		output_file
# 		]

# 		print(f"Running command: {' '.join(command)}")

# 		result = subprocess.run(command, capture_output=True, text=True)

# 		if result.returncode != 0:
# 			print(f"Error running SABREsim for {fn}")
# 			print("STDERR: ",result.stderr)
# 			continue

# 		energy = extractEfficienciesFromSTDOUT_kin4mc(fn)
# 		efficiencies = extractEfficienciesFromSTDOUT_kin4mc(result.stdout)

# 		if energy is not None:
# 			results.append({
# 					'energy':energy,
# 					'bu1':efficiencies['bu1'],
# 					'bu2':efficiencies['bu2'],
# 					'bu3':efficiencies['bu3'],
# 					'bu1_bu2':efficiencies['bu1_bu2'],
# 					'bu2_bu3':efficiencies['bu2_bu3'],
# 					'bu1_bu3':efficiencies['bu1_bu3'],
# 					'all':efficiencies['all']
# 				})
# 			count += 1

# 		print(f"bu1 eff = {efficiencies['bu1']:.4f}\nbu2 eff = {efficiencies['bu2']:.4f}\nbu3 eff = {efficiencies['bu3']:.4f}\nbu1_bu2 eff = {efficiencies['bu1_bu2']:.4f}\nbu2_bu3_eff = {efficiencies['bu2_bu3']:.4f}\nbu1_bu3 eff = {efficiencies['bu1_bu3']:.4f}\nall eff = {efficiencies['all']:.4f}\n")
# 		fn_counter += 1
# 		print(f"Finished {fn_counter} of {filenames_size}\t({fn_counter*100./filenames_size}%)")


# 	results.sort(key=lambda x:x['energy'])

# 	#output to a file! (I do NOT want to rerun SABREsim if I fuck up the plotting!)
# 	txt_output_path = os.path.join(output_dir, "kin4mc_efficiency_plot_data.txt")
# 	with open(txt_output_path, 'w') as f:
# 		f.write("ExE(keV)\tbu1\tbu2\tbu3\tbu1_bu2\tbu2_bu3\tbu1_bu3\tall")
# 		for i in range(len(energies)):
# 			f.write(f"{energies[i]}\t{bu1_vals[i]:.4f}\t{bu2_vals[i]:.4f}\t{bu3_vals[i]:.4f}\t{bu1_bu2_vals[i]:.4f}\t{bu2_bu3_vals[i]:.4f}\t{bu1_bu3_vals[i]:.4f}\t{all_vals[i]:.4f}")
# 	print(f"\n\nFinished running {count} files through SABREsim in {elapsed_time:.2f} seconds.")
# 	print(f"\n\nYou can find the data in {txt_output_path}\n")

def kin4mc(maxworkers=8):
	print("kin4mc() running...")
	start_time = time.time()

	input_dir = os.path.abspath(sys.argv[2])
	output_dir = os.path.abspath(sys.argv[3])
	os.makedirs(output_dir,exist_ok=True)

	filenames = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir,f))]
	filenames_size = len(filenames)

	results = []
	completed = 0

	print(f"Processing {filenames_size} files with {maxworkers} threads...")

	with ThreadPoolExecutor(max_workers=maxworkers) as executor:
		futures = {executor.submit(run_SABREsim_kin4mc, fn, input_dir, output_dir): fn for fn in filenames}

		for future in as_completed(futures):
			fn = futures[future]
			try:
				result = future.result()
				if result:
					results.append(result)
			except Exception as e:
				print(f"Error processing {fn}: {e}")
			finally:
				completed += 1
				print(f"Finished {completed}/{filenames_size} ({completed*100./filenames_size:.2f}%)")

	#sort by energy:
	results.sort(key=lambda x: x['energy'])

	#write to file:
	txt_output_path = os.path.join(output_dir,"kin4mc_efficiency_plot_data.txt")
	with open(txt_output_path, 'w') as f:
		f.write("ExE(keV)\tbu1\tbu2\tbu3\tbu1_bu2\tbu2_bu3\tbu1_bu3\tall\n")
		for r in results:
			f.write(f"{r['energy']:.4f}\t{r['bu1']:.4f}\t{r['bu2']:.4f}\t{r['bu3']:.4f}\t{r['bu1_bu2']:.4f}\t{r['bu2_bu3']:.4f}\t{r['bu1_bu3']:.4f}\t{r['all']:.4f}\n")

	elapsed_time = time.time() - start_time
	print(f"\nFinished running {len(results)} files through SABREsim in {elapsed_time:.2f} seconds!")
	print(f"Results saved to {txt_output_path}")


def main():
	#print("main function")
	numarg = len(sys.argv)
	if(numarg != 4):
		print("Invalid number of arguments! See calculateEfficiency.py for more info!")
		return

	if(sys.argv[1] == '2'):
		print("kin2mc chosen")
		print(f"path containing kin2mc.out files: {sys.argv[2]}")
		print(f"path containing SABREsim.det files: {sys.argv[3]}")
		kin2mc()
	elif(sys.argv[1] == '3'):
		print("kin3mc chosen")
		print(f"path containing kin3mc.out files: {sys.argv[2]}")
		print(f"path containing SABREsim.det files: {sys.argv[3]}")
		kin3mc()

	elif(sys.argv[1] == '4'):
		print("kin4mc")
		print(f"path containing kin4mc.out files: {sys.argv[2]}")
		print(f"path containing SABREsim.det files: {sys.argv[3]}")
		kin4mc()
	else:
		print("Invalid kinXmc argument! Expected 2,3,4 but got " + str(sys.argv[1]))
		return




if __name__=="__main__":
	main()
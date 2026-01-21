#this code takes the txt output from CalculateEfficiency and generates relative efficiency for all cases

import sys
import matplotlib.pyplot as plt


def kin2mc(filepath):
	#data storage:
	energies = []
	energymax = -1
	ej_vals = []
	ejmax = -1
	rec_vals = []
	recmax = -1
	both_vals = []
	bothmax = -1
	#first, let's read the file:
	with open(filepath, "r") as f:
		next(f)
		for line in f:
			parts = line.rstrip("\n").split("\t")

			energies.append(float(parts[0]))
			if(float(parts[0]) > energymax):
				energymax = float(parts[0])

			ej_vals.append(float(parts[1]))
			if(float(parts[1]) > ejmax):
				ejmax = float(parts[1])

			rec_vals.append(float(parts[2]))
			if(float(parts[2]) > recmax):
				recmax = float(parts[2])

			both_vals.append(float(parts[3]))
			if(float(parts[3]) > bothmax):
				bothmax = float(parts[3])

	ej_selfnorm = [e/ejmax for e in ej_vals]
	rec_selfnorm = [r/recmax for r in rec_vals]
	both_selfnorm = [b/bothmax for b in both_vals]

	plt.figure()
	plt.plot(energies, ej_selfnorm, label="Ejectile Only", marker=".")
	plt.xlabel("Excitation Energy (MeV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin2mc_relative_ej.png")
	plt.close()

	plt.figure()
	plt.plot(energies, rec_selfnorm, label="Recoil Only",marker=".")
	plt.xlabel("Excitation Energy (MeV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin2mc_relative_rec.png")
	plt.close()

	plt.figure()
	plt.plot(energies, both_selfnorm, label="Ejectile and Recoil", marker=".")
	plt.xlabel("Excitation Energy (MeV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin2mc_relative_both.png")
	plt.close()





def kin3mc(filepath):
	#data storage:
	energies = []
	energymax = -1
	ej_vals = []
	ejmax = -1
	bu1_vals = []
	bu1max = -1
	bu2_vals = []
	bu2max = -1
	both_vals = []
	bothmax = -1

	with open(filepath, "r") as f:
		next(f)
		for line in f:
			parts = line.rstrip("\n").split("\t")

			energies.append(float(parts[0]))
			if(float(parts[0]) > energymax):
				energymax = float(parts[0])

			bu1_vals.append(float(parts[1]))
			if(float(parts[1]) > bu1max):
				bu1max = float(parts[1])

			bu2_vals.append(float(parts[2]))
			if(float(parts[2]) > bu2max):
				bu2max = float(parts[2])

			both_vals.append(float(parts[3]))
			if(float(parts[3]) > bothmax):
				bothmax = float(parts[3])

	bu1_selfnorm = [b/bu1max for b in bu1_vals]
	bu2_selfnorm = [b/bu2max for b in bu2_vals]
	both_selfnorm = [b/bothmax for b in both_vals]

	plt.figure()
	plt.plot(energies, bu1_selfnorm, label="BU1 Only")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin3mc_relative_bu1.png")
	plt.close()

	plt.figure()
	plt.plot(energies, bu2_selfnorm, label="BU2 Only")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin3mc_relative_bu2.png")
	plt.close()

	plt.figure()
	plt.plot(energies, both_selfnorm, label="BU1 and BU2")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin3mc_relative_both.png")
	plt.close()

def kin4mc(filepath):
	#data storage:
	energies = []
	energymax = -1
	bu1_vals = []
	bu1max = -1
	bu2_vals = []
	bu2max = -1
	bu3_vals = []
	bu3max = -1
	bu1bu2_vals = []
	bu1bu2max = -1
	bu2bu3_vals = []
	bu2bu3max = -1
	bu1bu3_vals = []
	bu1bu3max = -1
	all_vals = []
	allmax = -1

	with open(filepath, "r") as f:
		next(f)
		for line in f:
			parts = line.rstrip("\n").split("\t")

			energies.append(float(parts[0]))
			if(float(parts[0]) > energymax):
				energymax = float(parts[0])

			bu1_vals.append(float(parts[1]))
			if(float(parts[1]) > bu1max):
				bu1max = float(parts[1])

			bu2_vals.append(float(parts[2]))
			if(float(parts[2]) > bu2max):
				bu2max = float(parts[2])

			bu3_vals.append(float(parts[3]))
			if(float(parts[3]) > bu3max):
				bu3max = float(parts[3])

			bu1bu2_vals.append(float(parts[4]))
			if(float(parts[4]) > bu1bu2max):
				bu1bu2max = float(parts[4])

			bu2bu3_vals.append(float(parts[5]))
			if(float(parts[5]) > bu2bu3max):
				bu2bu3max = float(parts[5])

			bu1bu3_vals.append(float(parts[5]))
			if(float(parts[5]) > bu1bu3max):
				bu1bu3max = float(parts[5])

			all_vals.append(float(parts[6]))
			if(float(parts[6]) > allmax):
				allmax = float(parts[6])

	bu1_selfnorm = [b/bu1max for b in bu1_vals]
	bu2_selfnorm = [b/bu2max for b in bu2_vals]
	bu3_selfnorm = [b/bu3max for b in bu3_vals]
	bu1bu2_selfnorm = [b/bu1bu2max for b in bu1bu2_vals]
	bu2bu3_selfnorm = [b/bu2bu3max for b in bu2bu3_vals]
	bu1bu3_selfnorm = [b/bu1bu3max for b in bu1bu3_vals]
	all_selfnorm = [a/allmax for a in all_vals]

	plt.figure()
	plt.plot(energies, bu1_selfnorm, label="BU1 Only")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin4mc_relative_bu1.png")
	plt.close()

	plt.figure()
	plt.plot(energies, bu2_selfnorm, label="BU2 Only")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin4mc_relative_bu2.png")
	plt.close()

	plt.figure()
	plt.plot(energies, bu3_selfnorm, label="BU3 Only")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin4mc_relative_bu3.png")
	plt.close()

	plt.figure()
	plt.plot(energies, bu1bu2_selfnorm, label="BU1 and BU2 Only")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin4mc_relative_bu1bu2.png")
	plt.close()

	plt.figure()
	plt.plot(energies, bu2bu3_selfnorm, label="BU2 and BU3 Only")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin4mc_relative_bu2bu3.png")
	plt.close()

	plt.figure()
	plt.plot(energies, bu1bu3_selfnorm, label="BU1 and BU3 Only")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin4mc_relative_bu1bu3.png")
	plt.close()

	plt.figure()
	plt.plot(energies, all_selfnorm, label="BU1, BU2, BU3 Only")
	plt.xlabel("Excitation Energy (keV)")
	plt.ylabel("Relative Efficiency (%)")
	plt.ylim(0,1)
	plt.title("Relative Detection Efficiency vs Excitation Energy")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.savefig("/mnt/e/SABREsim/efficiencies/kin4mc_relative_all.png")
	plt.close()

def main():
	numarg = len(sys.argv)
	if(numarg != 3):
		print("Invalid number of arguments! See CalculateEfficiency.py for more info!")
		return

	if(sys.argv[1] == '2'):
		print("kin2mc option chosen")
		kin2mc(sys.argv[2])

	elif(sys.argv[1] == '3'):
		print("kin3mc option chosen")
		kin3mc(sys.argv[2])

	elif(sys.argv[1] == '4'):
		print("kin4mc option chosen")
		kin4mc(sys.argv[2])


if __name__ == "__main__":
	main()
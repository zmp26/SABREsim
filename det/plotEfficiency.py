import pandas as pd
import matplotlib.pyplot as plt




def plotKin4mc():
	filename = "/home/zmpur/efficiencies/kin4mc/Boron9_p8Be_4He4He/dets/kin4mc_efficiency_plot_data.txt"

	df = pd.read_csv(filename, delim_whitespace=True, comment='#', engine='python')
	print(df.head())

	plt.figure(figsize=(10,6))

	# plt.plot(df["ExE(keV)"], df["bu1"], label="bu1")
	# plt.plot(df["ExE(keV)"], df["bu2"], label="bu2")
	# plt.plot(df["ExE(keV)"], df["bu3"], label="bu3")
	# plt.plot(df["ExE(keV)"], df["bu1_bu2"], label="bu1_bu2")
	# plt.plot(df["ExE(keV)"], df["bu2_bu3"], label="bu2_bu3")
	# plt.plot(df["ExE(keV)"], df["bu1_bu3"], label="bu1_bu3")
	# plt.plot(df["ExE(keV)"], df["all"], label="all")\

	df_norm = df.copy()
	df_norm["bu1"] = df["bu1"]/df["bu1"].max()
	df_norm["bu2"] = df["bu2"]/df["bu2"].max()
	df_norm["bu3"] = df["bu3"]/df["bu3"].max()
	df_norm["bu1_bu2"] = df["bu1_bu2"]/df["bu1_bu2"].max()
	df_norm["bu2_bu3"] = df["bu2_bu3"]/df["bu2_bu3"].max()
	df_norm["bu1_bu3"] = df["bu1_bu3"]/df["bu1_bu3"].max()
	df_norm["all"] = df["all"]/df["all"].max()

	plt.plot(df["ExE(keV)"], df_norm["all"], label="all")


	plt.xlabel("ExE (keV)")
	plt.ylabel("Efficiency")
	plt.title("Efficiencies for 10B(3He,4He)9B, 9B->8Be+p, 8Be->4He+4He")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.show()

def plotKin3mc():
	filename = "/home/zmpur/efficiencies/kin3mc/Lithium6_ad/dets/kin3mc_efficiency_plot_data.txt"

	df = pd.read_csv(filename, delim_whitespace=True, comment='#', engine='python')
	print(df.head())

	plt.figure(figsize=(10,6))

	# plt.plot(df["ExE(keV)"], df["bu1"], label="bu1")
	# plt.plot(df["ExE(keV)"], df["bu2"], label="bu2")
	# plt.plot(df["ExE(keV)"], df["both"], label="both")

	df_norm = df.copy()
	df_norm["bu1"] = df["bu1"]/df["bu1"].max()
	df_norm["bu2"] = df["bu2"]/df["bu2"].max()
	df_norm["both"] = df["both"]/df["both"].max()

	plt.plot(df["ExE(keV)"], df_norm["both"], label="both")

	plt.xlabel("ExE (keV)")
	plt.ylabel("Efficiency")
	plt.title("Efficiencies for 7Li(3He,4He)6Li, 6Li->a+d")
	plt.legend()
	plt.grid(True)
	plt.tight_layout()
	plt.show()
output_file = "temp.out"

bu1_eff = 5.932
bu2_eff = 4.905
both_eff = 88.849

with open(output_file, "w") as f:
	for energy in range(0, 1500, 5):
		f.write(f"{energy:.1f}\t{bu1_eff}\t{bu2_eff}\t{both_eff}\n")

print(f"Data written to {output_file}")
import shutil
import re

filename = "/mnt/e/TRANSMIT_6LI_3061keV_LiF_4000A.txt"
backup_filename = filename + ".bak"

shutil.copy(filename,backup_filename)

with open(filename, "r") as f:
	lines = f.readlines()


newlines = []
for line in lines:
	# if line and line[0] in "TBS" and len(line) > 1 and line[1].isdigit():
	# 	line = line[0] + " " + line[1:]
	# 	if line[1] == 0 and line[2] == 0:
	if line and line[0] in "TBS":
		m = re.match(r"([TBS])(\d+)(.*)",line)
		if m:
			type_char = m.group(1)
			ion_number = m.group(2)
			rest = m.group(3)
			line = f"{type_char} {ion_number}{rest}\n"

	newlines.append(line)


with open(filename,"w") as f:
	f.writelines(newlines)

print(f"File {filename} has had its spacing fixed!")
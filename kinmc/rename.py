import os
import re

directory = "/mnt/e/SABREsim/kinmc/kin4mc_eff_10B3He4He9B_p_4He"
'''
kin4mc_7500keV_10B3He4He9B_at0_keV_p_4He_4He.out --> kin4mc_7500keV_10B3He4He9B_at0keV_p_4He.out
'''

for filename in os.listdir(directory):
    if filename.endswith(".out"):
        filenamesplit = filename.split("_")
        newfilename = filenamesplit[0] + "_" + filenamesplit[1] + "_" + filenamesplit[2] + "_" + filenamesplit[3] + "keV_" + filenamesplit[5] + "_" + filenamesplit[6] + ".out"
        #print(f"old filename = {filename}\tnew filename = {newfilename}\n")
        old_path = os.path.join(directory,filename)
        new_path = os.path.join(directory,newfilename)

        os.rename(old_path,new_path)
        print(f"Renaming {filename} to {newfilename}")

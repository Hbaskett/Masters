# script to read the NASA9 coefficients from vulcan and write them into a file for use in ATMO

from pathlib import Path
import numpy as np

NASA9_dir = Path.home()/"vulcan"/"thermo"/"NASA9"
ATMO_NASA9 = Path.home()/"atmo"/"chem"

# these species have 6 rows of coefficients
rogues_1 = ["C2H4OH_1","C2H4OH_2","CH2CHO","HNO3","CH3O2","CH3OCO","C4H9_2","C2H3CHO","C2H5O","SH","CO2H","HS","C2H5OO","C2H4HO2","CH3ONO","HOCN","HCNO","NH2OH","CH3NO2","HO2","HONO","C2H5OOH","HNOH","HONO2","C4H9_1","C2H6CO","H2CNO","C3H7_1","H2NO","HCNH","C2H5CHO","NaOH","HCNN"]
# this species has 2 rows of coefficients
rogues_2 = ["C4H9"]

ATMO_names = ["cC6H6","O-3P","O-1D","CHCO","SN","OOH","1CH2","N-4S","3CH2"]
VULCAN_names = ["C6H6","O","O_1","HCCO","NS","HO2","CH2_1","N","CH2"]

#**/*.txt for all .txt in NASA9_dir

for p in sorted(Path(NASA9_dir).glob('**/*.txt')):
    #print(f"{p.name}:\n{p.read_text()}\n")
    
    # get the name of the file
    fname = p.name
    species = fname.removesuffix(".txt")
    
    # some species have completely different naming convernsions so rename them
    for i in range(0,len(VULCAN_names)):
        if species == VULCAN_names[i]:
            species = species.replace(species,ATMO_names[i])
    
    # replace _ with - to format to python
    for letter in species:
        if letter == "_":
            species = species.replace(letter,"-")
    
    # for the standard format of 4 lines of polynomials
    with open(NASA9_dir/fname) as f:
        vul_lines = f.readlines()
            
    f.close()
    
    line = ["",""]
    
    if species not in (rogues_1 and rogues_2):
    
        # split each element of vul_lines into induvidual strings
        coeffs_1 = vul_lines[0].split()
        coeffs_2 = vul_lines[1].split()
        coeffs_3 = vul_lines[2].split()
        coeffs_4 = vul_lines[3].split()
        line=["",""]
    
        # write each string + space into line array as a single line
        for i in range(0,len(coeffs_1)):
            if i == 0:
                if coeffs_1[0].startswith("-"):
                    line[0] += coeffs_1[i]
                else:
                    line[0] += " " + coeffs_1[i]
            
                if coeffs_3[0].startswith("-"):
                    line[1] += coeffs_3[i]
                else:
                    line[1] += " " + coeffs_3[i]    
        
            if i != 0:
                if coeffs_1[i].startswith("-"):
                    line[0] += " " + coeffs_1[i]
                else:
                    line[0] += "  " + coeffs_1[i]
            
                if coeffs_3[i].startswith("-"):
                    line[1] += " " + coeffs_3[i]
                else:
                    line[1] += "  " + coeffs_3[i]
                
        for i in range(0,len(coeffs_2)):
            # omit coefficients that are 0
            if i != 2:
                if coeffs_2[i].startswith("-"):
                    line[0] += " " + coeffs_2[i]
                else:
                    line[0] += "  " + coeffs_2[i]
            
                if coeffs_4[i].startswith("-"):
                    line[1] += " " + coeffs_4[i]
                else:
                    line[1] += "  " + coeffs_4[i]
    
    # formating not uniform in the NASA9 VULCAN files so introduce if statements to deal with rogues
    
    # this species has 2 rows of coefficients
    if species == "C4H9":
    
        coeffs_1 = vul_lines[0].split()
        coeffs_2 = vul_lines[1].split()
        line=["",""]
        
        for i in range(0,len(coeffs_1)):
            if i == 0:
                if coeffs_1[0].startswith("-"):
                    line[0] += coeffs_1[i]
                else:
                    line[0] += " " + coeffs_1[i]
            
                if coeffs_2[0].startswith("-"):
                    line[1] += coeffs_2[i]
                else:
                    line[1] += " " + coeffs_2[i]
                
            if i != 0:
                # omit coefficients that are 0
                if i != 2:
                    if coeffs_1[i].startswith("-"):
                        line[0] += " " + coeffs_1[i]
                    else:
                        line[0] += "  " + coeffs_1[i]
                
                    if coeffs_2[i].startswith("-"):
                        line[1] += " " + coeffs_2[i]
                    else:
                        line[1] += "  " + coeffs_2[i]
    
    # these species have 6 rows of coefficients - program chooses the bottom 4 rows
    if species in rogues_1:
    
        #coeffs_1 = vul_lines[0].split()
        #coeffs_2 = vul_lines[1].split()
        coeffs_3 = vul_lines[2].split()
        coeffs_4 = vul_lines[3].split()
        coeffs_5 = vul_lines[4].split()
        coeffs_6 = vul_lines[5].split()
        
        line=["",""]
        
        for i in range(0,len(coeffs_3)):
            if i == 0:
                if coeffs_3[0].startswith("-"):
                    line[0] += coeffs_3[i]
                else:
                    line[0] += " " + coeffs_3[i]
            
                if coeffs_5[0].startswith("-"):
                    line[1] += coeffs_5[0]
                else:
                    line[1] += " " + coeffs_5[0]    
        
            if i != 0:
                if coeffs_3[i].startswith("-"):
                    line[0] += " " + coeffs_3[i]
                else:
                    line[0] += "  " + coeffs_3[i]
            
                if coeffs_5[i].startswith("-"):
                    line[1] += " " + coeffs_5[i]
                else:
                    line[1] += "  " + coeffs_5[i]
                
        for i in range(0,len(coeffs_4)):
            # omit coefficients that are 0
            if i != 2:
                if coeffs_4[i].startswith("-"):
                    line[0] += " " + coeffs_4[i]
                else:
                    line[0] += "  " + coeffs_4[i]
            
                if coeffs_6[i].startswith("-"):
                    line[1] += " " + coeffs_6[i]
                else:
                    line[1] += "  " + coeffs_6[i]
      
    with open(ATMO_NASA9/"coeff9_NASA_sc_vulcan.dat", mode = "a") as f:
    
        #need to formate the number 200 to the 22nd column of the .dat file so make case by case for molecule lengths
        if len(species) == 1: # 20 spaces
            f.write(f"{species}                    200  6000  1000\n")
        if len(species) == 2: # 19 spaces
            f.write(f"{species}                   200  6000  1000\n")
        if len(species) == 3:
            f.write(f"{species}                  200  6000  1000\n")
        if len(species) == 4:
            f.write(f"{species}                 200  6000  1000\n")
        if len(species) == 5:
            f.write(f"{species}                200  6000  1000\n")
        if len(species) == 6:
            f.write(f"{species}               200  6000  1000\n")
        if len(species) == 7:
            f.write(f"{species}              200  6000  1000\n")
        if len(species) == 8:
            f.write(f"{species}             200  6000  1000\n")
        if len(species) == 9:
            f.write(f"{species}            200  6000  1000\n")
        if len(species) == 10:
            f.write(f"{species}           200  6000  1000\n")
        if len(species) == 11:
            f.write(f"{species}          200  6000  1000\n")
        if len(species) == 12:
            f.write(f"{species}         200  6000  1000\n")
        if len(species) == 13:
            f.write(f"{species}        200  6000  1000\n")
        if len(species) == 14:
            f.write(f"{species}       200  6000  1000\n")
        f.write(f"{line[1]}\n")
        f.write(f"{line[0]}\n")
        
    f.close()
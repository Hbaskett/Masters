# script to convert Lodders abundances to ATMO format
# list of A_x in form from Lodders 2009 found in vulcan/fastchem_vulcan/input/solar_elemental_abundance.dat
# A_x = np.log(n_x/n_H)+12 abundance from Lodders 2009
A = {"C" : 8.4434,
      "H" : 12.00,
      "He" : 10.9864,
      "N" : 7.9130,
      "O" : 8.7826,
      "P" : 5.5058,
      "S" : 7.12,
      "Si" : 7.5867,
      "Ti" : 4.9794,
      "V" : 4.0437,
      "Cl" : 5.3002,
      "K" : 5.1619,
      "Na" : 6.3479,
      "Mg" : 7.5995,
      "F" : 4.49196,
      "Ca" : 6.3677,
      "Fe" : 7.5151,
}


# normalised to value of silicon
n_si = 1E6

# then calculate n_x/n_H for ATMO
H_ratio = {}
for ele in A:
    x = 10**(A[ele]-12)
    H_ratio = {**H_ratio, ele: x}
    
print(H_ratio)

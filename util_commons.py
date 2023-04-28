# Python file containing all the useful dictionaries e.t.c. for scripts
from pathlib import Path
from collections import OrderedDict

# dictionary of linestlyes
linestyles_dict = OrderedDict(
    [('solid',               (0, ())),
    #('loosely dotted',      (0, (1, 10))),
    ('dotted',              (0, (1, 5))),
    #('densely dotted',      (0, (1, 1))),

    #('loosely dashed',      (0, (5, 10))),
    ('dashed',              (0, (5, 5))),
    #('densely dashed',      (0, (5, 1))),

    #('loosely dashdotted',  (0, (3, 10, 1, 10))),
    ('dashdotted',          (0, (3, 5, 1, 5))),
    #('densely dashdotted',  (0, (3, 1, 1, 1))),

    ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
    #('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
    ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

simulations = { "atmo" : {
                    "fdir" : Path.home()/"atmo"/"mphys_final_runs",
                    "trans_dir" : Path.home()/"atmo"/"mphys_final_runs"/"trans_runs",
                    "exp" : [
                    #"HD189_EQ_atherm_ncon", # ATMO NASA9 coefficients 
                    #"HD189_EQ_vtherm_ncon", # VULCAN NASA9 coefficients
                    #"HD189_EQ_inversion_ncon",
                    #"HD189_default_kzz1e9_ncon", # only kzz and mol
                    #"HD189_default_kzz1e9_inversion_ncon",
                    "HD189_full_kzz1e9_ncon", # vz = 0
                    #"HD189_full_kzz1e9_inversion_ncon",
                    "HD189_full_kzz1e9_vz1e1_ncon",
                    "HD189_full_kzz1e9_vz1e2_ncon",
                    "HD189_full_kzz1e9_vz1e3_ncon",
                    "HD189_full_kzz1e9_vz1e4_ncon",
                    "HD189_full_kzz1e9_vz1e5_ncon",
                    #"HD189_full_kzz1e9_vz1e6_ncon",
                    #"HD189_full_kzz1e9_vz1e7_ncon",
                    #"HD189_full_kzz1e9_vz1e8_ncon",
                    "HD189_full_kzz1e9_vz-1e1_ncon",
                    "HD189_full_kzz1e9_vz-1e2_ncon",
                    "HD189_full_kzz1e9_vz-1e3_ncon",
                    "HD189_full_kzz1e9_vz-1e4_ncon",
                    "HD189_full_kzz1e9_vz-1e5_ncon",
                    #"HD189_full_kzz1e9_vz-1e6_ncon",
                    #"HD189_full_kzz1e9_vz-1e7_ncon",
                    #"HD189_full_kzz1e9_vz-1e8_ncon"
                    ]
                    },
                    
                "vulcan" : {
                    "fdir" : Path.home()/"vulcan"/"output",
                    "exp" : [
                    #"HD189_EQ", # turned off all forcing in the _cfg file
                    #"HD189_default_kzz1e9", #only kzz and mol
                    #"HD189_full_kzz1e9", # vz = 0
                    #"HD189_full_kzz1e9_inversion",
                    #"HD189_full_kzz1e9_vz1e1",
                    #"HD189_full_kzz1e9_vz1e2",
                    #"HD189_full_kzz1e9_vz1e3",
                    #"HD189_full_kzz1e9_vz1e4",
                    #"HD189_full_kzz1e9_vz1e5",
                    #"HD189_full_kzz1e9_vz1e6",
                    #"HD189_full_kzz1e9_vz1e7",
                    #"HD189_full_kzz1e9_vz1e8",
                    #"HD189_full_kzz1e9_vz-1e1",
                    #"HD189_full_kzz1e9_vz-1e2",
                    #"HD189_full_kzz1e9_vz-1e3",
                    #"HD189_full_kzz1e9_vz-1e4",
                    #"HD189_full_kzz1e9_vz-1e5",
                    #"HD189_full_kzz1e9_vz-1e6",
                    #"HD189_full_kzz1e9_vz-1e7",
                    #"HD189_full_kzz1e9_vz-1e8"
                    ]
                 },            
    }
    
# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180),(255, 127, 14),(44, 160, 44),(214, 39, 40),(148, 103, 189),(140, 86, 75), (227, 119, 194),(127, 127, 127),(188, 189, 34),(23, 190, 207),(174, 199, 232),(255, 187, 120),(152, 223, 138),(255, 152, 150),(197, 176, 213),(196, 156, 148),(247, 182, 210),(199, 199, 199),(219, 219, 141),(158, 218, 229)]

# Species list for chem_tsai2020.ncdf
species_list = ['H2','H','O-3P','OH','H2O','CH','C','3CH2','CH3','CH4','C2','C2H2','C2H3','C2H','C2H4','C2H5','C2H6','CO','CO2','CH2OH','H2CO','HCO','CH3O','CH3OH','CH3CO','O2','CH2CO','CHCO','He','N-4S','NH','CN','HCN','NO','NH2','N2','NH3','N2H2','N2H','N2H3','N2H4','HNO','H2CN','HNCO','NO2','N2O','C4H2','CH2NH2','CH2NH','CH3NH2','CH3CHO','HNO2','NCO','OOH','H2O2','HC3N','CH3CN','CH2CN','C2H3CN','SH','HSO','H2S','C3H3', 'C3H2', 'C3H4', 'C6H5', 'cC6H6','S','S2','SO','CS','COS','CS2','SN','HS2','SO2','S4','S8','HCS','S3','C4H3','C4H5','S2O','CH3SH','CH3S','O_1D','1CH2','N_2D']
    
# Species list from chem_funs.py in vulcan
spec_list = ['OH','H2','H2O','H','O','CH','C','CH2','CH3','CH4','C2','C2H2','C2H','C2H3','C2H4','C2H5','C2H6','CO','CO2','CH2OH','H2CO','HCO','CH3O','CH3OH','CH3CO','O2','H2CCO','HCCO','N','NH','CN','HCN','NO','NH2','N2','NH3','N2H2','N2H','N2H3','N2H4','HNO','H2CN','HNCO','NO2','N2O','C4H2','CH2NH2','CH2NH','CH3NH2','CH3CHO','HNO2','NCO','HO2','H2O2','HC3N','CH3CN','CH2CN','C2H3CN','SH','HSO','H2S','C3H3','C3H2','C3H4','C6H5','C6H6','S','S2','SO','CS','COS','CS2','NS','HS2','SO2','S4','S8','HCS','S3','C4H3','C4H5','S2O','CH3SH','CH3S','O_1','CH2_1','N_2D','He']
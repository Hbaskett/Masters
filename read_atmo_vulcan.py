# Script which plots atmo and vulcan abundances on the same set of axis
import iris
import warnings
from netCDF4 import *
from constant import *
from pathlib import Path
from pylab             import *
from numpy             import *
from math              import *
from scipy.interpolate import *
from scipy.signal      import *
from scipy.optimize    import *
from scipy.special     import sph_harm
from scipy.integrate   import *
from aeolus.plot import add_custom_legend, subplot_label_generator
import matplotlib              as mpl
import matplotlib.font_manager as mpfm
import mpl_toolkits.mplot3d as m3d
from collections import OrderedDict
from tqdm.notebook import tqdm as tqdm


try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False
import os, sys
import pickle
from util_commons import tableau20, simulations, spec_list, species_list, linestyles_dict

warnings.filterwarnings("ignore", module="iris")
warnings.filterwarnings("ignore", module="aeolus")

mpl.rc('text',usetex=False)    
mpl.rc('font',family='serif',size=18)

prop = mpfm.FontProperties(size=12)

def plot_atmo_vulcan(atmo_imol=[1,18,5,10,2,29,37], figname=""):

    # Parameters
    for_legend_simu = {} #dictionary of simulations for labels
    color_index=0 # colour index for clr dictionary
    clr = {} # dictionary of colours
    color_index_legend = 0 # colour index for legend
    for_legend_gas = {} # dictionary of species for labels
    for_simu = {} # contains linestyles
   
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
    for i in range(len(tableau20)):    
        r, g, b = tableau20[i]    
        tableau20[i] = (r / 255., g / 255., b / 255.)
    
    # resize atmo_imol
    if atmo_imol==[]:
        atmo_imol = array(range(nmol))
    else:
        atmo_imol = array(atmo_imol)-1
    
    vrbls={}
    
    for model in ["atmo","vulcan"]:
        
        vrbls[model] = {}
        fdir = simulations[model]["fdir"]
        
        for simu in simulations[model]["exp"]:
    
            vrbls[model][simu] = {}
            if model == "atmo":
            
                atmo_file = simu + '.ncdf'
                file = fdir/atmo_file

                ncdfFile = Dataset(file,'r')
                vars = ncdfFile.variables

                nlevel = len(ncdfFile.dimensions['nlevel'])
                ab = vars['abundances'][:,:]
                pp = vars['pressure'][:]/1.e6
                nm = vars['molname'][:,:]
                
                # produce array of labels for the species which will be plotted
                species = []
                for i in range(len(atmo_imol)):
                    molname = nm[atmo_imol[i],:]  
                    molname = molname.tostring() # convert nm to a string
                    molname = molname.decode() # remove 'b' prefix from string
                    molname = molname.strip() # remove extra spaces after letters in string
                    species = np.append(species,molname)
            
                    # write pressure and abundance to dictionary
                    vrbls[model][simu][molname] = {
                            "pressure": pp,
                            "abundances": ab[atmo_imol[i],:],
                    }
                    
            if model == "vulcan":

                vulcan_file = simu + '.vul'
                vul_file = fdir/vulcan_file
        
                with open(vul_file, 'rb') as handle:
                  data = pickle.load(handle)
        
                ppr = data['atm']['pco']/1.e6
                molecules = data['variable']['species']

                # write pressure and abundance to dictionary
                for i in range(len(atmo_imol)):
                    molname = species_list[atmo_imol[i]]  #find the molecule name from the Tsai chemical list
              
                    if simu == "HD189_EQ":
                        vrbls[model][simu][molname] = {
                            "pressure": data['atm']['pco']/1.e6,
                            "abundances": data['variable']['ymix'][:,molecules.index(molname)][::1],
                        }
                    else:
                        vrbls[model][simu][molname] = {
                            "pressure": data['atm']['pco']/1.e6,
                            "abundances": data['variable']['ymix'][:,molecules.index(molname)],
                        }
    
    for gas in species:
        clr[gas] = {"color": tableau20[color_index]}
        color_index = color_index + 1
        for_legend_gas.update({gas:{"linewidth": 1.5,"color": tableau20[color_index_legend]}})
        color_index_legend =  color_index_legend + 1
    
    # update dictionaries with linestyles for plotting
    iterable = iter(linestyles_dict)
    
    for model in ["atmo","vulcan"]:
        for simu in simulations[model]["exp"]:
            x = next(iterable)
        
            # make array to store run name parameters
            run = []
            #run.append(model)
        
            # split simu into parts
            legend_name = simu.split('_')
        
            # extract the planet name from the file name
            planet = legend_name[0]
        
            # extract the type of run
            #if "full" in simu:
            #    run.append(legend_name[1])
        
            if "default" in simu:
                run.append(legend_name[1])
        
            if "EQ" in simu:
                run.append(legend_name[1])
                # NASA9 coefficients
                #run.append(legend_name[2])
        
            # extract the kzz value
            if "kzz" in simu:
                kzz_value = legend_name[2]
                kzz_power = int(kzz_value[-1])
                # write kzz value to legend
                #if kzz_value[3] == "-":
                #    run.append(f"k_zz = -{10**kzz_power}") 
                #else:
                #    run.append(f"k_zz = {10**kzz_power}")
        
            # extract the vz value
            if simu == ("HD189_full_kzz1e9_ncon"):
                run.append("$v_z$ = 0" " [$cm$$s^{-1}$]")
            
            # extract the vz value
            if simu == ("HD189_full_kzz1e9"):
                run.append("$v_z$ = 0 [$cm$$s^{-1}$]")
            
            if "vz" in simu:
                vz_value = legend_name[3]
                vz_power = int(vz_value[-1])
                # write vz value to legend
                if vz_value[2] == "-":
                    run.append(f"$v_z$ = -{10**vz_power}" " [$cm$$s^{-1}$]")
                else:
                    run.append(f"$v_z$ = {10**vz_power}" " [$cm$$s^{-1}$]")
        
            if "inversion" in simu:
                run.append("inversion")
        
            full = ""
            for element in run:
                full += element + " "
        
            for_simu.update({simu: {"linestyle": linestyles_dict[x]}})
            for_legend_simu.update({full: {"linestyle": linestyles_dict[x],"linewidth": 1, "color": "k"}}) 
    
    # plot species for each simulation
    fig, ax = plt.subplots(
        ncols=1, nrows=1, figsize=(14,9), sharex=True, sharey=True, constrained_layout=True,
    )
    iletters = subplot_label_generator()

    for model in ["atmo","vulcan"]:
        for simu in simulations[model]["exp"]:
            for gas in species:
                ax.plot(
                    vrbls[model][simu][gas]["abundances"].data,
                    vrbls[model][simu][gas]["pressure"].data,
                    **clr[gas],
                    **for_simu[simu],
                    linewidth=1.5,
                )
                
    # Add extra legends
    #add_custom_legend(
    #    ax, for_legend_simu,bbox_to_anchor=(0,0.40),loc="center left", frameon=False, fontsize=12
    #)
    
    add_custom_legend(
        ax, for_legend_gas,bbox_to_anchor=(0,0),loc="lower left", frameon=False, fontsize=15
    )
               
    # Remove duplicate labels
    handles, labels = plt.gca().get_legend_handles_labels()
    i = 1
    while i < len(labels):
        if labels[i] in labels[:i]:
            del labels[i]
            del handles[i]
        else:
            i += 1  
                
    # Apply legend 
    leg = ax.legend(
        handles,
        labels,
        loc="upper right",
        bbox_to_anchor=(-0.04, 1),
        frameon=False,
        handlelength=0,
        fontsize="small",
    )
            
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid()
    ylim(1e3,1e-9) 
    xlim(1E-12,1E0)  
    ax.set_xlabel('Abundances', fontsize=20)
    ax.set_ylabel('Pressure [bar]', fontsize=20)
    plt.savefig(figname + ".png",bbox_inches="tight", dpi=600)
    
plot_atmo_vulcan(atmo_imol=[18,5,10,2,37], figname="tester")

# 1 H2
# 29 He
# 18 CO
# 5 H2O
# 10 CH4
# 2 H
# 37 NH3
# Script which plots atmo and vulcan transport flux on the same set of axis
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

#import vulcan_cfg
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
mpl.rc('font',family='serif',size=12)
mpl.rc('xtick', labelsize=10)
mpl.rc('ytick', labelsize=10)

prop = mpfm.FontProperties(size=12)

def plot_phi(atmo_imol=[1,4,9,12,34,62,68]):
     
    # Parameters
    ns = 88 # number of species in network
    diff_comps = ["diff_mol","diff_kzz","diff_buoy","diff_therm","diff_vz"] # names of diffusion terms 
    vrbls = {} # dictionary of data for plotting
    species = [] # array of species for plotting
    color_index=0 # colour index for clr dictionary
    clr = {} # dictionary of colours for simualtions
    for_simu = {}
    for_legend_simu = {} # dictionary of simulations labels
    for_legend_gas = {} # dictionary of species labels
    for_legend_phi = {} # dictionary for phi labels
          
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
    for i in range(len(tableau20)):    
        r, g, b = tableau20[i]    
        tableau20[i] = (r / 255., g / 255., b / 255.)
    
    # resize atmo_imol
    if atmo_imol==[]:
        atmo_imol = array(range(nmol))
    else:
        atmo_imol = array(atmo_imol)-1
        
    for model in ["atmo","vulcan"]:
    
        vrbls[model] = {}
        fdir = simulations[model]["fdir"]
        
        for simu in simulations[model]["exp"]:
        
            vrbls[model][simu] = {}
            if model == "atmo":
                
                chem_fname = simu + ".ncdf"
                phi_fname = simu + "_atmo_phi_comps" + ".dat"
                file = fdir/chem_fname
    
                ncdfFile = Dataset(file,'r')
                vars = ncdfFile.variables
                
                nlevel = 151
                pp = vars['pressure'][:]
                nm = vars['molname'][:,:]
                
                #conversion in bar
                pp = pp/1.0E6
                
                # load phi components into array from data file
                phi_comps = np.genfromtxt(fdir/phi_fname,comments = "=",usecols = (0,1,2,3,4))
                
                # produce array of labels for the species which will be plotted
                for i in range(len(atmo_imol)):
                    molname = nm[atmo_imol[i],:]  
                    molname = molname.tostring() # convert nm to a string
                    molname = molname.decode() # remove 'b' prefix from string
                    molname = molname.strip() # remove extra spaces after letters in string
                    species = np.append(species,molname)
                    #print(molname)
                    
                    # write pressure to plotting dictionary
                    vrbls[model][simu][molname] = {
                            "pressure": pp[0:nlevel],
                    }
        
                    diff_total = np.zeros(shape=(nlevel)) # set total transport flux to 0
                    
                    for j in range(len(diff_comps)):
                        diff_plot = np.zeros(shape=(ns,nlevel)) # define as a column of zeros
                        diff_plot = phi_comps[:,j] # slice the jth component of phi into diff_phi
                    
                        # for some reason when vz = 0 get negligable values of adv_comp which mess up plotting so default it to 0
                        if (j == 4 and simu == "HD189_full_kzz1e9_ncon"):
                            diff_plot = np.zeros(shape=(ns*nlevel))
                    
                        diff_plot = diff_plot.reshape(ns,nlevel) # reshape into usable shape
                        diff_plot = diff_plot[atmo_imol[i]:atmo_imol[i]+1,:] # take slice for the ith species
                        diff_plot = diff_plot.reshape(nlevel)
                        diff_total += diff_plot # update total phi component
                        diff_plot = abs(diff_plot) # take absolute values
                        diff_total = abs(diff_total)
                        
                        # write to plotting dictionary
                        vrbls[model][simu][molname].update({diff_comps[j] : diff_plot})
                        vrbls[model][simu][molname].update({"diff_total" : diff_total})
            
            if model == "vulcan":
                # find the simulation file
                chem_fname_vulcan = simu + '.vul'
                chem_file_vulcan = fdir/chem_fname_vulcan
                # find the simulation phi components file
                phi_fname_vulcan = simu + "_phi_comps" + ".dat"
                phi_file_vulcan = fdir/phi_fname_vulcan
                
                # open the simulation file and extact the pressure data and list of molecules
                with open(chem_file_vulcan, 'rb') as handle:
                  data = pickle.load(handle)
                  
                ppr = data['atm']['pco']/1.e6
                molecules = data['variable']['species']
                nz = 150 # number of vertical layers from vulcan_cfg
        
                # open the phi components file
                with open(phi_file_vulcan, 'rb') as handle: 
                  phi_data = pickle.load(handle)
                
                # make an array of the names of the species that will be plotted
                for i in range(len(atmo_imol)):
                      molname = species_list[atmo_imol[i]]  #find the molecule name from the Tsai chemical list
                      #print(molname)
                      y = spec_list.index(molname) # select the index of the species from the spec list
                      
                      # write pressure to plotting dictionary
                      vrbls[model][simu][molname] = {
                          "pressure": ppr,
                      }
                
                      diff_total = np.zeros(nz) # set total phi to 0 for all levels
                      for j in diff_comps:
                      
                          diff_plot = np.zeros(nz) # set plotting diff to 0 for all levels
                          diff_comp = phi_data[j] # slice the jth component 
                          #if j != "diff_vz":
                          #    diff_comp = phi_data[j] * 1e6 # artificial units scaling m^2 to cm^2
                          diff_plot = diff_comp[:,y:y+1].flatten() # flattern shape of diff_comp array to a column
                          diff_total += diff_plot # add to total phi
                          diff_plot = abs(diff_plot) # take absolute values
                          diff_total = abs(diff_total)
                          
                          # write phi components to plotting dictionary 
                          vrbls[model][simu][molname].update({j : diff_plot})
                          vrbls[model][simu][molname].update({"diff_total" : diff_total})
                          
    # update dictionaries with linestyles for plotting
    iterable = iter(linestyles_dict) # for multiple species
    linestyles = ["solid","dotted","dashed"] # for 1 species
    
    for model in ["atmo","vulcan"]:
        
        #i = 0 # use for 1 species
        for simu in simulations[model]["exp"]:
            
            #x = next(iterable) # use for multiple gases
            run = [] # make array to store run name parameters
        
            # Append the model the simulation was run on
            run.append(model)
        
            # split simu into parts
            legend_name = simu.split('_')
        
            # extract the planet name from the file name
            planet = legend_name[0]
        
            # extract the type of run 
            run.append(legend_name[1])
        
            # extract the vz value
            if simu == "HD189_full_kzz1e9_ncon":
                run.append("$v_z$ = 0 [$cm$$s^{-1}$]")
        
            if simu == "HD189_full_kzz1e9":
                run.append("$v_z$ = 0 [$cm$$s^{-1}$]")
        
            if "vz" in simu:
                vz_value = legend_name[3]
                vz_power = int(vz_value[-1])
            
                if vz_value[2] == "-":
                    run.append(f"$v_z$ = -{10**vz_power}" " [$cm$$s^{-1}$]")
                else:
                    run.append(f"$v_z$ = {10**vz_power}" " [$cm$$s^{-1}$]")
        
            full = ""
            for element in run:
                full += element + " "
            
            #for_simu.update({simu: {"linestyle": linestyles_dict[x]}}) # for multiple species
            #for_legend_simu.update({full: {"linestyle": linestyles_dict[x],"linewidth": 1, "color": "k"}})
            
            for_simu.update({simu: {"linestyle": "solid", "color": tableau20[color_index]}}) # for 1 species
            for_legend_simu.update({full: {"linestyle": linestyles[i],"linewidth": 1.5, "color": tableau20[color_index]}})
            color_index = color_index + 1
            #i = i + 1
    
    # use when plotting multiple gases
    # add colors for each species to dictionary for legend
    for gas in species:
        clr[gas] = {"color": tableau20[color_index]}
        #for_legend_gas.update({gas:{"linewidth": 1.5,"color": "k"}}) # for 1 species
        for_legend_gas.update({gas:{"linewidth": 1.5,"color": tableau20[color_index]}})
        #color_index = color_index + 1
    
    diff_comps.append("diff_total") # append diff total to diff comps to make plotting easier
    
    # Plot
    fig, axes = plt.subplots(
        ncols=2, nrows=1, figsize=(12,6), sharex=True, sharey=True, constrained_layout=True,
    )
    iletters = subplot_label_generator()
    
    #colour = "r"
    for model in ["atmo","vulcan"]:
    
        #if model == "atmo":
        #    colour = "r"
        #if model == "vulcan":
        #    colour = "k"
    
        for simu in simulations[model]["exp"]:
            for gas in species:
                for comp, ax in zip(["diff_vz","diff_total"], axes.flatten()):
                    for_legend_phi = {comp: {"linestyle": linestyles_dict["solid"],"linewidth": 1, "color": "k"}}
                    ax.plot(
                        vrbls[model][simu][gas][comp].data,
                        vrbls[model][simu][gas]["pressure"].data,
                        #color = colour,
                        #**clr[gas], # use for multiple species
                        **for_simu[simu],
                        linewidth=1.5,
                    )
                    
                    if comp == "diff_total":
                        # add gas legend in top right plot
                        add_custom_legend(
                            ax, for_legend_gas, loc="upper right", frameon=False, fontsize=12
                        )
                        
                        # Add simulations legend
                        #add_custom_legend(
                        #    ax, for_legend_simu, bbox_to_anchor=(0.7,0.8), loc="upper right", frameon=False, fontsize=12
                        #) 
                    
                    # Add phi components legend
                    #add_custom_legend(
                    #    ax, for_legend_phi, bbox_to_anchor=(1,0.90), loc="upper right", frameon=False, fontsize=20
                    #)
    
                    ax.set_yscale('log')
                    ax.set_xscale('log')
                    ax.grid(visible="true")
                    ylim(1e3,1e-9)
                    xlim(1E-6,1E20) 
                    # for only 5 plots
                    ax.set_xlabel('Magnitude of Transport Flux [$cm^{-2}$$s^{-1}$]', fontsize=12)
                    ax.set_ylabel('Pressure [bar]', fontsize=12)
                     
    #plt.setp(axes[-1, :], xlabel='Magnitude of Transport Flux [$cm^{-2}$$s^{-1}$]')
    #plt.setp(axes[:, 0], ylabel='Pressure [bar]')
    
    # to plot only 5 subplots
    #axes[2][1].set_visible(False)
    #axes[2][0].set_position([0.25,0,0.5,0.31])
              
    plt.savefig(f"test",bbox_inches="tight", dpi=300)
    
    
plot_phi(atmo_imol=[2])

# 1 H2
# 29 He
# 18 CO
# 5 H2O
# 10 CO
# 2 H
# 37 NH3
from pathlib import Path
from pylab import *
from netCDF4 import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker
import iris
import warnings
from netCDF4 import *
from constant import *
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
from util_commons import tableau20, simulations, spec_list, species_list, linestyles_dict

def plot_trans(figname=""):


    # Parameters
    for_legend_simu = {} # dictionaty of simualtions for labels
    for_simu = {} # dictionary for simulations
    model = "atmo"
    color_index = 0
    
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
    for i in range(len(tableau20)):    
        r, g, b = tableau20[i]    
        tableau20[i] = (r / 255., g / 255., b / 255.)
    
    vrbls = {}
    if model == "atmo":
    
        vrbls[model] = {}
        for simu in simulations[model]["exp"]:
        
            atmo_file = simu + '_trans' + '.ncdf'
            file = simulations[model]["trans_dir"]/atmo_file
                
            spectrum = Dataset(file)
    
            nu = spectrum.variables['nu'][:]
            rprs = spectrum.variables['transit_radius'][:]
            
            #Convert to 'transit depth'
            #rprs = rprs**2.
            
            nu = 1./nu
            nu = nu/100.
            nu = nu*1E6
            
            vrbls[model][simu] = {
                      "wavelength" : nu[:],
                      "radius ratio" : rprs[:],
            }
    
    # update dictionaries with linestyles for plotting
    iterable = iter(linestyles_dict)
    
    for model in ["atmo","vulcan"]:
        for simu in simulations[model]["exp"]:
            x = next(iterable)
        
            # make array to store run name parameters
            run = []
            run.append(model)
        
            # split simu into parts
            legend_name = simu.split('_')
        
            # extract the planet name from the file name
            planet = legend_name[0]
        
            # extract the type of run
            if "full" in simu:
                run.append(legend_name[1])
        
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
        
            for_simu.update({simu: {"linestyle": "solid", "color": tableau20[color_index]}})
            for_legend_simu.update({full: {"linestyle": "solid","linewidth": 1.5, "color": tableau20[color_index]}})
            color_index = color_index + 1
    
    # plot species for each simulation
    fig, ax = plt.subplots(
        ncols=1, nrows=1, figsize=(14,9), sharex=True, sharey=True, constrained_layout=True,
    )
    iletters = subplot_label_generator()

    for model in ["atmo","vulcan"]:
        for simu in simulations[model]["exp"]:
            semilogx(
                vrbls[model][simu]["wavelength"].data,
                vrbls[model][simu]["radius ratio"].data,
                **for_simu[simu],
                linewidth=1.5,
            )
                        
    #Add extra legends
    add_custom_legend(
        ax, for_legend_simu,bbox_to_anchor=(0,0.80),loc="center left", frameon=False, fontsize=12
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
         
    xlim(0.3,30)  
    ax.set_xlabel('Wavelength $\mu$m', fontsize=20)
    ax.set_ylabel('R$_p$/R$_s$', fontsize=20)
    
    #Set xaxis ticker format
    ax = gca()
    majorLocator = FixedLocator([0.4,0.6,1.0,2,4,8,16,24])
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    tick_params(axis='x',labelsize='20')
    tick_params(axis='y',labelsize='20')
    
    plt.savefig(figname + ".png",bbox_inches="tight", dpi=600)

plot_trans(figname="test")
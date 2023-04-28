# Python script to plot PT profiles from ATMO

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

import sys
#sys.path.insert(0, '../') # including the upper level of directory for the path of modules

try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False
import os, sys
import pickle

warnings.filterwarnings("ignore", module="iris")
warnings.filterwarnings("ignore", module="aeolus")

mpl.rc('text',usetex=False)    
mpl.rc('font',family='serif',size=12)

mpl.rc('xtick', labelsize=20)
mpl.rc('ytick', labelsize=20)

prop = mpfm.FontProperties(size=12)

def plot_pt(name=""):

    # dictionary of linestyles
    linestyles_dict = OrderedDict(
        [('solid',               (0, ())),
         ('loosely dotted',      (0, (1, 10))),
         ('dotted',              (0, (1, 5))),
         ('densely dotted',      (0, (1, 1))),

         ('loosely dashed',      (0, (5, 10))),
         ('dashed',              (0, (5, 5))),
         ('densely dashed',      (0, (5, 1))),

         ('loosely dashdotted',  (0, (3, 10, 1, 10))),
         ('dashdotted',          (0, (3, 5, 1, 5))),
         ('densely dashdotted',  (0, (3, 1, 1, 1))),

         ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
         ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
         ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
         
    # Set up colors for each species   
    # These are the "Tableau 20" colors as RGB.    
    tableau20 = [(31, 119, 180),(255, 127, 14),(44, 160, 44),(214, 39, 40),(148, 103, 189),(140, 86, 75), (227, 119, 194),(127, 127, 127),(188, 189, 34),(23, 190, 207),(174, 199, 232),(255, 187, 120),(152, 223, 138),(255, 152, 150),(197, 176, 213),(196, 156, 148),(247, 182, 210),(199, 199, 199),(219, 219, 141),(158, 218, 229)] 
    
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
    for i in range(len(tableau20)):    
        r, g, b = tableau20[i]    
        tableau20[i] = (r / 255., g / 255., b / 255.)
     
    pt_names = {"HD189_PT_vulcan_isotherm_extended" : "a" ,
                "HD189_PT_vulcan_inversion_extended" : "a" ,
    }

    vrbls={}
    for simu in pt_names:

        fdir = Path.home()/"atmo"/"mphys_final_runs"
        fname = simu + ".ncdf"
    
        file = fdir/fname

        ncdfFile = Dataset(file,'r')
        vars = ncdfFile.variables

        tt = vars['temperature'][:]
        pp = vars['pressure'][:]

        #conversion in bar
        pp = pp/1.0E6

        vrbls[simu] = {
                "pressure" : pp,
                "temperature" : tt,
        }


    color_index_legend = 0
    for_legend_pt = {}
    for simu in pt_names:
        
        run = []
        # split simu into parts
        legend_name = simu.split('_')
        
        #planet = legend_name[0]
        
        version = legend_name[3]
        run.append(version)

        full = ""
        for element in run:
            full += element + " "
            
        # add colours for each species to dictionary for legend
        for_legend_pt.update({full:{"linewidth": 1,"color": tableau20[color_index_legend]}})
        color_index_legend =  color_index_legend + 1
    
    # plot species for each simulation
    fig, ax = plt.subplots(
        ncols=1, nrows=1, figsize=(14,9), sharex=True, sharey=True, constrained_layout=True,
    )
    iletters = subplot_label_generator()

    for simu, line in pt_names.items():
        ax.plot(
            vrbls[simu]["temperature"].data,
            vrbls[simu]["pressure"].data,
            linewidth=1,
        )
    
    #add_custom_legend(
    #    ax, for_legend_pt,bbox_to_anchor=(0,0),loc="lower left", frameon=False, fontsize=12
    #)
    
    ax.set_yscale('log')    
    #ax.set_title("HD189733b PT profiles",fontsize=20)
    ylim((1e3,1e-9))
    ax.set_xlabel('Temperature [K]', fontsize=20)
    ax.set_ylabel('Pressure [bar]', fontsize=20)
    plt.savefig(name + ".png",bbox_inches="tight", dpi=600)

plot_pt(name="HD189_PT_profiles")
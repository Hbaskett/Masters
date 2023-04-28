from pylab import *
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
import matplotlib              as mpl
import matplotlib.font_manager as mpfm
import mpl_toolkits.mplot3d as m3d
import netCDF4                 as nc

def extend_pt(fp="",fn="",name="",nfig=1,tit='',clear=False,label='',cond = '',color = '',linewidth = 1.,linestyle='-',marker='',new_ncdf=""):

    fdir = Path.home()/"atmo"
    fpath = fdir/fp
    fname = fn + ".ncdf"
    figname = name + ".png"
    
    file = fpath/fname
    
    # read pressure and temperature from the PT ncdf file
    ncdfFile = Dataset(file,'r')
    vars = ncdfFile.variables

    tt = vars['temperature'][:]
    pp = vars['pressure'][:]
    
    ncdfFile.close()
    
    #interpolate new temperatures in the upper atmosphere up to 10-8 bar by adding an isotherm at the same temperature of the final temperature value and just extend the atmosphere upwards
    
    # Extract the final value of temperature at the top of the atmosphere
    tt_final = tt[0]
    
    #reverse the arrays so the upper atmosphere values are at the "end"
    tt_rev = np.flip(tt)
    pp_rev = np.flip(pp)
    
    #append the isotherm into the array
    tt_extend_rev = np.append(tt_rev,tt_final)
    pp_extend_rev = np.append(pp_rev,1E-3)
    
    #reverse the array back to there original setup
    tt_extend = np.flip(tt_extend_rev)
    pp_extend = np.flip(pp_extend_rev)
    
    # Write new ncdf file with extended pressures and temperatures
    new_ncdf = new_ncdf + ".ncdf"
    
    new_ncdf = nc.Dataset(new_ncdf,'w',format='NETCDF3_CLASSIC')

    ndepth = pp_extend.size

    if ((not pp_extend.size == ndepth) or (not tt_extend.size == ndepth)):
        print('the arrays must have the same size')
        return

    new_ncdf.createDimension('nlevel',ndepth)
    ntt=new_ncdf.createVariable('temperature','f8',('nlevel',))
    npp=new_ncdf.createVariable('pressure','f8',('nlevel',))
    nkzz=new_ncdf.createVariable('kzz','f8',('nlevel',))
    nvz=new_ncdf.createVariable('vz','f8',('nlevel',))
    
    #if not kzz == 'none':
    #    nkk=nout.createVariable('kzz','f8',('nlevel',))
    #    nkk.units = 'cm2 s-1'
    #    nkk[:] = kzz[:]

    ntt.units  = 'K'
    npp.units  = 'dyn cm-2'
    nkzz.units = 'cm2 s-1'
    nvz.units = 'cm s-1'
    
    kzz = np.linspace(1000000000000,33200000,num=ndepth)
    vz = np.linspace(10000,5000,num=ndepth)
    
    npp[:]  = pp_extend[:]
    ntt[:]  = tt_extend[:]
    nkzz[:] = kzz[:]
    nvz[:] = vz[:]

    new_ncdf.close()
    
    #Plot the extended PT structure
    #conversion in bar
    pp = pp/1.0E6

    figure(nfig,figsize=(15,10))
    if clear:
        clf()
    
    if (cond == 'h2o' or cond == 'H2O' or cond == 'both' or cond == 'BOTH') :
    	    semilogy(10000./(38.84-3.83*numpy.log10(pp_extend)-3.93*0.-0.20*0.*numpy.log10(pp_extend)),pp_extend,label = 'H2O condensation line',color='red',lw=linewidth)
    if (cond == 'nh3' or cond == 'NH3' or cond == 'both' or cond == 'BOTH') :
    	    semilogy(10000./(68.02-6.31*numpy.log10(pp_extend)-6.19*0.),pp_extend,label = 'NH3 condensation line',color='magenta',lw=linewidth)
   	    
    if (color=='') :
    	    semilogy(tt_extend, pp_extend,label=label,lw=linewidth,linestyle=linestyle,marker=marker)
    else :
    	    semilogy(tt_extend, pp_extend,label=label,color=color,lw=linewidth,linestyle=linestyle,marker=marker)
    	    
    title(tit,fontsize=12)
    ylim((pp_extend.max(),pp_extend.min()))
    xlabel('Temperature [K]')
    ylabel('Pressure [bar]')

    if not label=='':
        legend(prop=prop).draw_frame(0)
    show()
    
extend_pt(fp="mphys_final_runs",fn="HD189_PT_vulcan_inversion",name="HD189_PT_vulcan_inversion_extended",new_ncdf="HD189_PT_vulcan_inversion_extended")
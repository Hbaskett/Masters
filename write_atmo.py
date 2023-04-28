from pylab             import *
from numpy             import *
from math              import *
import matplotlib              as mpl
import matplotlib.font_manager as mpfm
import netCDF4                 as nc

mpl.rc('text',usetex=False)    
mpl.rc('font',family='serif',size=12)

prop = mpfm.FontProperties(size=12)

def write_atmo(pp,tt,kzz='none',fout='input.ncdf'):

    nout = nc.Dataset(fout,'w',format='NETCDF3_CLASSIC')

    ndepth = pp.size

    if ((not pp.size == ndepth) or (not tt.size == ndepth)):
        print 'the arrays must have the same size'
        return

    nout.createDimension('nlevel',ndepth)
    ntt=nout.createVariable('temperature','f8',('nlevel',))
    npp=nout.createVariable('pressure','f8',('nlevel',))
    if not kzz == 'none':
        nkk=nout.createVariable('kzz','f8',('nlevel',))
        nkk.units = 'cm2 s-1'
        nkk[:] = kzz[:]

    ntt.units  = 'K'
    npp.units  = 'dyn cm-2'
    

    npp[:]  = pp[:]
    ntt[:]  = tt[:]

    nout.close()

def write_abundances(fin='inidat_dummy',fout='chem_dummy.ncdf',nlevel=64,pmin=1E-5,pmax=1E3):

    #dummy pressure

    pp=zeros(nlevel)
    pp[0]      = pmin*1E6 #bar -> cgs
    pp[nlevel-1] = pmax*1E6 #bar -> cgs

    dlp = (log10(pmax)-log10(pmin))/(nlevel-1)

    for i in range(1,nlevel,1):
        pp[i] = pp[i-1]*10.**dlp

    #elemental abundances David 1996
    aH  = 0.909964  
    aHe = 0.088714  
    aC  = 3.26E-04  
    aN  = 1.02E-04  
    aO  = 4.77E-04  
    aTi = 0.#8.13E-08  
    aV  = 0.#7.76E-09  

    nH2  = aH/2.
    nHe  = aHe
    nCH4 = aC
    nN2  = aN/2.
    nO2  = aO/2.
    nTiO = aTi
    nVO  = aV

    ntot = nH2+nHe+nCH4+nN2+nO2+nTiO+nVO

    nH2  = nH2 /ntot 
    nHe  = nHe /ntot
    nCH4 = nCH4/ntot
    nN2  = nN2 /ntot
    nO2  = nO2 /ntot
    nTiO = nTiO/ntot
    nVO  = nVO /ntot  

    nin = loadtxt(fin,dtype=np.string_)

    molname = nin[:,0]
    nmol = nin[:,0].size
    lname = 10

    Amol = zeros((nlevel,nmol))

    nH2 = 0.853
    nHe = 0.145
    nN2 = 7.11e-5
    nCH4 = 5.66e-4
    nO2 = 5.78e-4

    for i in range(nmol):
        if nin[i,0]=='H2':
            Amol[:,i] = nH2
        elif nin[i,0] == 'He':
            Amol[:,i] = nHe
        elif nin[i,0] == 'CH4':
            Amol[:,i] = nCH4
        elif nin[i,0] == 'N2':
            Amol[:,i] = nN2
        elif nin[i,0] == 'O2':
            Amol[:,i] = nO2
        elif nin[i,0] == 'TiO':
            Amol[:,i] = nTiO
        elif nin[i,0] == 'VO':
            Amol[:,i] = nVO

   # for i in range(nlevel):
   #     Amol[i,:] = nin[:,1]

    
    molname_blanck = array([' '])
    for i in range(lname-1):
        molname_blanck = np.append(molname_blanck,' ')

    molname_line = np.copy(molname_blanck)
    molname_line[:len(nin[0,0])] = list(nin[0,0])
    molname = array([molname_line])
    
    for j in range(1,nmol):
        molname_line = np.copy(molname_blanck)
        molname_line[:len(nin[j,0])] = list(nin[j,0])
        molname = np.append(molname,[molname_line],axis=0)

    molname = transpose(molname)

    nout = nc.Dataset(fout,'w',format='NETCDF3_CLASSIC')
    nout.createDimension('nlevel',nlevel)
    nout.createDimension('lname',lname)
    nout.createDimension('nmol',nmol)

    #dummy pressure 
    npp=nout.createVariable('pressure','f8',('nlevel',))
    npp[:] = pp[:]

    # mol name
    nmolname=nout.createVariable('molname','c',('nmol','lname',))
    nmolname[:,:] = transpose(molname[:,:] )
    

    # Abundances
    namol=nout.createVariable('abundances','f8',('nmol','nlevel',))
    namol[:,:] = transpose(Amol[:,:] )

    

    nout.close()
    

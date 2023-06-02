"""
Collection of codes to calculate stellar flux collected by the 
istrument.

@author: A. G. Sreejith
"""

#########################################
###    Import Libraries and Functions
#########################################
import os
import numpy as np
import astropy.modeling.functional_models as am
import csc_functions as csc 
from astropy.io import ascii
import scipy.special as ss
import scipy.constants as sc
import astropy.constants as ac

#Constants
    R_sun = 6.957e10                     # cm
    AU    = 1.49598e13                   # cm (Astronomical Unit)
    c0    = 2.99792458e10                # cm/s (speed of light)
    sigma = 5.67051e-5                   # erg/cm**2/s/K**4 (stefan-boltzmann)
    k_B   = 1.380658e-16                 # erg/K = g cm**2/s**2/K (boltzmann const)
    N_A   = 6.02214179e23                # /mol (Avagadro constant)
    cA    = 2.99792458e18                # Angstrom/s (speed of light)
    vc    = ac.c.to("km/s").value        # Speed of light, km/s
    #===============================

class Voigt(object):
    """
    1D Voigt profile model.

    Parameters
    ----------
    x0: Float
       Line center location.
    hwhmL: Float
       Half-width at half maximum of the Lorentz distribution.
    hwhmG: Float
       Half-width at half maximum of the Gaussian distribution.
    scale: Float
       Scale of the profile (scale=1 returns a profile with integral=1.0).

    Example
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.opacity.broadening as b
    >>> Nl = 5
    >>> Nw = 10.0
    >>> hG = 1.0
    >>> HL = np.logspace(-2, 2, Nl)
    >>> l = b.Lorentz(x0=0.0)
    >>> d = b.Gauss  (x0=0.0, hwhm=hG)
    >>> v = b.Voigt  (x0=0.0, hwhmG=hG)

    >>> plt.figure(11, (6,6))
    >>> plt.clf()
    >>> plt.subplots_adjust(0.15, 0.1, 0.95, 0.95, wspace=0, hspace=0)
    >>> for i in np.arange(Nl):
    >>>   hL = HL[i]
    >>>   ax = plt.subplot(Nl, 1, 1+i)
    >>>   v.hwhmL = hL
    >>>   l.hwhm  = hL
    >>>   width = 0.5346*hL + np.sqrt(0.2166*hL**2+hG**2)
    >>>   x = np.arange(-Nw*width, Nw*width, width/1000.0)
    >>>   plt.plot(x/width, l(x), lw=2.0, color="b",         label="Lorentz")
    >>>   plt.plot(x/width, d(x), lw=2.0, color="limegreen", label="Doppler")
    >>>   plt.plot(x/width, v(x), lw=2.0, color="orange",    label="Voigt",
    >>>            dashes=(8,2))
    >>>   plt.ylim(np.amin([l(x), v(x)]), 3*np.amax([l(x), v(x), d(x)]))
    >>>   ax.set_yscale("log")
    >>>   plt.text(0.025, 0.75, r"$\rm HW_L/HW_G={:4g}$".format(hL/hG),
    >>>            transform=ax.transAxes)
    >>>   plt.xlim(-Nw, Nw)
    >>>   plt.xlabel(r"$\rm x/HW_V$", fontsize=12)
    >>>   plt.ylabel(r"$\rm Profile$")
    >>>   if i != Nl-1:
    >>>       ax.set_xticklabels([""])
    >>>   if i == 0:
    >>>       plt.legend(loc="upper right", fontsize=11)
    """
    def __init__(self, x0=0.0, hwhmL=1.0, hwhmG=1.0, scale=1.0):
        # Profile parameters:
        self.x0    = x0
        self.hwhmL = hwhmL
        self.hwhmG = hwhmG
        self.scale = scale
        # Constants:
        self._A = np.array([-1.2150, -1.3509, -1.2150, -1.3509])
        self._B = np.array([ 1.2359,  0.3786, -1.2359, -0.3786])
        self._C = np.array([-0.3085,  0.5906, -0.3085,  0.5906])
        self._D = np.array([ 0.0210, -1.1858, -0.0210,  1.1858])
        self._sqrtln2 = np.sqrt(np.log(2.0))
        self._sqrtpi  = np.sqrt(np.pi)


    def __call__(self, x):
        return self.eval(x)


    def eval(self, x):
        """
        Compute Voigt profile over the specified coordinates range.

        Parameters
        ----------
        x: 1D float ndarray
           Input coordinates where to evaluate the profile.

        Returns
        -------
        v: 1D float ndarray
           The line profile at the x locations.
        """
        if self.hwhmL/self.hwhmG < 0.1:
            sigma = self.hwhmG / (self._sqrtln2 * np.sqrt(2))
            z = (x + 1j * self.hwhmL - self.x0) / (sigma * np.sqrt(2))
            return self.scale * ss.wofz(z).real / (sigma * np.sqrt(2*np.pi))

        # This is faster than the previous script (but fails for HWl/HWg > 1.0):
        X = (x-self.x0) * self._sqrtln2 / self.hwhmG
        Y = self.hwhmL * self._sqrtln2 / self.hwhmG

        V = 0.0
        for i in np.arange(4):
            V += (self._C[i]*(Y-self._A[i]) + self._D[i]*(X-self._B[i])) \
                 / ((Y-self._A[i])**2 + (X-self._B[i])**2)
        V /= np.pi * self.hwhmL
        return \
            self.scale * self.hwhmL/self.hwhmG * self._sqrtpi*self._sqrtln2 * V




def voigtq(wavelength, absorber, line):
    bnorm = absorber["B"] / vc
    # Doppler width (Hz):
    vd = absorber["B"]*sc.kilo / (line["wave"] * sc.angstrom)

    vel = (wavelength/(line["wave"]*(1.0+absorber["Z"])) - 1.0)  / bnorm
    a = line["gamma"] / (4*np.pi * vd)

    idx_wings = np.abs(vel) >= 10.0
    idx_core = np.abs(vel) < 10.0

    vo = vel*0.0
    if np.any(idx_wings) > 0:
          vel2 = vel[idx_wings]**2
          hh1 = 0.56419/vel2 + 0.846/vel2**2
          hh3 = -0.56/vel2**2
          vo[idx_wings] = a * (hh1 + a**2 * hh3)

    if np.any(idx_core) > 0:
        x0 = 0.0
        hwhm_L = line["gamma"] / np.sqrt(2) / 2
        hwhm_G = vd * np.sqrt(np.log(2))
        voigt = Voigt(x0, hwhm_L, hwhm_G)
        vo[idx_core] = voigt(vel[idx_core]*vd)
        vo[idx_core] /= np.amax(vo[idx_core])

    tau = 0.014971475*(10.0**absorber["N"]) * line["F"] * vo/vd

    return np.exp(-tau)


def gaussian(wavelength, wl0, sigma, scale=1.0):
  """
  Compute a Gaussian function centered at wl0, with width sigma, and
  height scale at wl0.
  """
  return scale * np.exp(-0.5*((wavelength-wl0)/sigma)**2)

def cute_snr_lca(flux,sigmaMg22,sigmaMg21,t_star,r_star,logr,stype):
    
    # Scaling relation for stellar types
    #  F5 V-G9 V 16  Mg II 14.1-267     -0.291  -0.0208 24.2  31.2
    #  K0 V-K5 V 12  Mg II 8.4-86.6     0.338   -0.318  29.8  39.3
    #  M0 V-M5 V 07  Mg II 0.032-17.6   0.814   -0.296  27.8  33.8
    #  F5 V-G9 V 13  Ca II 11.4-149.    -0.223  0.080   42.2  49.6
    #  K0 V-K5 V 10  Ca II 5.9-39.4     0.605   -0.433  20.6  30.7
    #  M0 V-M5 V 07  Ca II 0.0076-6.53  1.028   -0.312  39.9  49.6

    #spectral values for interpolation
#    old data
#    R_book = [0.58,0.61,0.65,0.68,0.7,0.73,0.75,0.78,0.81,0.84,0.88,0.91,0.93,0.96,0.97,0.99,
#              1.02,1.04,1.08,1.11,1.16,1.21,1.26,1.32,1.36]
#    T_book = [3935,4045,4258,4557,4791,4925,5047,5156,5273,5386,5484,5560,5626,5678,5723,
#              5767,5819,5870,5948,6017,6115,6226,6332,6445,6540]   # T_eff of book
#    BV_book = [1.48,1.42,1.31,1.15,1.03,0.966,0.912,0.866,0.819,0.776,0.740,0.713,0.69,0.672,
#               0.657,0.642,0.625,0.608,0.583,0.561,0.530,0.496,0.464,0.431,0.404]

    data=ascii.read(os.path.join(os.path.curdir,'extra', 'stellar_param_mamjeck.txt'),comment='#')
        
    Sp=data['col1']
    T_book=data['col2']
    R_book=data['col9']
    BV_book=data['col8']
    Teff=float(t_star)
    loc=csc.find_nearest(T_book,Teff)
   
    try:
        stype
    except NameError:
        stype  = Sp[loc] 
    try:
        BV
    except NameError:
        BV  = BV_book[loc] 
    try:
        r_star
    except NameError:
        r_star  = R_book[loc]

    logR=logr
    if (stype == 'F5V' or stype == 'F6V' or stype == 'F7V' or 
        stype == 'F8V' or stype == 'F9V' or stype == 'F9.5V' or 
        stype == 'G0V' or stype == 'G1V' or stype == 'G2V' or 
        stype == 'G3V' or stype == 'G4V' or stype == 'G5V' or
        stype == 'G6V' or stype == 'G7V' or stype == 'G8V' or 
        stype == 'G9V') :
                c1 = 0.87
                c2 = 5.73
                Rmg = 10**(c1*logR+c2)
    elif (stype == 'K9V' or stype == 'K8V' or stype == 'K7V' or 
          stype == 'K6.5V' or  stype == 'K6V' or stype == 'K5.5V' or 
          stype == 'K5V' or stype == 'K4.5V' or stype == 'K4V' or 
          stype == 'K3.5V' or stype == 'K3V' or stype == 'K2.5V' or
          stype == 'K2V' or stype == 'K1.5V' or stype == 'K1V' or 
          stype == 'K0.5V' or stype == 'K0V') :
                c1 = 1.01
                c2 = 6.00
                Rmg = 10**(c1*logR+c2)
    elif (stype == 'M5V' or stype == 'M4V' or stype == 'M3V' or 
          stype == 'M2V' or stype == 'M2.5V' or stype == 'M1.5V' or 
          stype == 'M1V' or stype == 'M0.5V' or stype == 'M0V')  :
                c1 = 1.59
                c2 = 6.96
                Rmg = 10**(c1*logR+c2)
    else:
                Rmg = 0.0
    #===============================
    #Parameters for MgII
    
    MgII1w      = 2795.5280 #MgIIh wavelength 
    MgII1_loggf = 0.100     #MgIIh loggf
    MgII1_stark = -5.680    #MgIIh stark
    MgII2w      = 2802.7050 #MgIIk wavelength
    MgII2_loggf = -0.210    #MgIIk loggf
    MgII2_stark =-5.680     #MgIIkStark damping constant

    Mgaratio_loggf2to1=(10**MgII2_loggf)/(10**MgII1_loggf)
    ##===============================
    in_flux=flux[:,1]/flux[:,2]
    fl=flux[:,1]
    ##===============================
    #Parameters for CaII
    #TCa2=7.7d5             #temperature of the gas at the formation of the Ca2 emission in K
    #mCa2=40.078/N_A        #atomic mass of Ca2 in g: atomic mass in g/mole / Avogrado number
    #CaHwl0=3968.470        #CaH wavelength
    #CaH_loggf=-0.2         #CaH loggf
    #CaH_stark=8.193        #Stark damping constant
    #
    #CaKwl0=3933.664        #CaK wavelength
    #CaK_loggf=0.105        #CaK loggf
    #CaK_stark=8.207        #Stark damping constant
    #
    #Caratio_loggfKtoH=10**CaK_loggf/10**CaH_loggf
    ##ISM fixed parameters
    #ISM_b_Mg2=2.0          #b-parameter for the Ca2 ISM lines in km/s
    #vr_ISM=0.              #radial velocity of the ISM absorption lines in km/s
    ##===============================
    
    E=Rmg*AU**2
    #E=50*AU**2
    Mg21em=E/(1.+Mgaratio_loggf2to1)
    Mg22em=Mgaratio_loggf2to1*E/(1.+Mgaratio_loggf2to1)
    
    gaussMg22=gaussian(flux[:,0],MgII2w,sigmaMg22,0.3989*Mg22em/sigmaMg22)
    gaussMg21=gaussian(flux[:,0],MgII1w,sigmaMg21,0.3989*Mg21em/sigmaMg21)
    '''
    CaII calculations
    #Ca2em=E/(1.+Caratio_loggfKtoH)
    #Ca2Kem=Caratio_loggfKtoH*E/(1.+Caratio_loggfKtoH)
    #
    #sigmaCa2K=sqrt(k_B*TCa2/mCa2/c0**2)*CaKwl0# & sigmaCa2K=8.*sigmaCa2K
    #gaussCa2K=gaussian(flux[0,*],[0.3989*Ca2Kem/sigmaCa2K,CaKwl0,sigmaCa2K])
    #
    #sigmaCa2H=sqrt(k_B*TCa2/mCa2/c0**2)*CaHwl0# & sigmaCa2H=8.*sigmaCa2H
    #gaussCa2H=gaussian(flux[0,*],[0.3989*Ca2Hem/sigmaCa2H,CaHwl0,sigmaCa2H])
    '''
    gaussMg2 = gaussMg21 + gaussMg22
    flux_emission = flux[:,1] + gaussMg2
    flux[:,1]= flux_emission
    return flux
    

def cute_ism_abs_all(flux,n_mg2,n_mg1,n_fe2):
    #Constants
    R_sun = 6.957e10       # cm
    AU    = 1.49598e13     # cm (Astronomical Unit)
    vc    = 299792.458e0   #speed of light in km/s
    c0    = 2.99792458e10  # cm/s (speed of light)
    sigma = 5.67051e-5     # erg/cm**2/s/K**4 (stefan-boltzmann)
    k_B   = 1.380658e-16   # erg/K = g cm**2/s**2/K (boltzmann const)
    N_A   = 6.02214179e23  # /mol (Avagadro constant)
    cA    = 2.99792458e18  # Angstrom/s (speed of light)
    #===============================

    ##Parameters for MgII
    MgII1w      = 2795.5280
    MgII1_loggf = 0.100
    MgII1_stark = -5.680
    MgII2w      = 2802.7050
    MgII2_loggf = -0.210
    MgII2_stark =-5.680     
    Mgaratio_loggf2to1=(10**MgII2_loggf)/(10**MgII1_loggf)
  
 
    ##Parameters for MgI
    MgIw      = 2852.127 
    MgI_loggf = 0.255
    MgI_stark = -5.640
    
    ##Parameters for FeII
    FeIIw      = 2599.39515
    FeII_loggf = 0.378
    FeII_stark = -6.53
    
    #ISM fixed parameters
    ISM_b_Mg2=2.0        #b-parameter for the Ca2 ISM lines in km/s
    vr_ISM=0.            #radial velocity of the ISM absorption lines in km/s
    

    #Construct the ISM absorption
    n_flux=flux[:,1]/flux[:,2]
    
    #for MgII doublet
    absorberMg1={'ion':'MG21','N': n_mg2,'B':ISM_b_Mg2,'Z': 0.0}
    lineMg1={'ion':'Mg21','wave':MgII1w+MgII1w*vr_ISM/vc,'F':10**MgII1_loggf,'gamma':10**MgII1_stark}
    ISMMg21=voigtq(flux[:,0],absorberMg1,lineMg1)
    
    absorberMg2={'ion':'MG22','N':n_mg2,'B':ISM_b_Mg2,'Z':0.0}
    lineMg2={'ion':'Mg22','wave':MgII2w+MgII2w*vr_ISM/vc,'F':10**MgII2_loggf,'gamma':10**MgII2_stark}
    ISMMg22=voigtq(flux[:,0],absorberMg2,lineMg2)
    
    #for MgI
    absorberMgI={'ion':'MG1','N':n_mg1,'B':ISM_b_Mg2,'Z':0.0}
    lineMgI={'ion':'Mg1','wave':MgIw+MgIw*vr_ISM/vc,'F':10**MgI_loggf,'gamma':10**MgI_stark}
    ISMMg1=voigtq(flux[:,0],absorberMgI,lineMgI)
    
    #for FeII 
    absorberFeII={'ion':'FE2"','N':n_fe2,'B':ISM_b_Mg2,'Z':0.0}
    lineFeII={'ion':'Fe2','wave':FeIIw+FeIIw*vr_ISM/vc,'F':10**FeII_loggf,'gamma':10**FeII_stark}
    ISMFe2=voigtq(flux[:,0],absorberFeII,lineFeII)
    
    #for CaII
    #absorberCaK=create_struct('ion','Ca2K','N',N,'B',ISM_b_Ca2,'Z',0.0)
    #lineCaK=create_struct('ion','Ca2K','wave',CaKwl0+CaKwl0*vr_ISM/vc,'F',10**CaK_loggf,'gamma',10**CaK_stark)
    #ISMCaK=voigtq(flux[0,*],absorberCaK,lineCaK)
    #
    #absorberCaH=create_struct('ion','Ca2H','N',N,'B',ISM_b_Ca2,'Z',0.0)
    #lineCaH=create_struct('ion','Ca2H','wave',CaHwl0+CaHwl0*vr_ISM/vc,'F',10**CaH_loggf,'gamma',10**CaH_stark)
    #ISMCaH=voigtq(flux[0,*],absorberCaH,lineCaH)
    
    ISM = ISMMg21*ISMMg22*ISMMg1*ISMFe2
    
    flux_absorption = ISM * n_flux
    flux[:,1]=flux[:,2]*flux_absorption
    return flux

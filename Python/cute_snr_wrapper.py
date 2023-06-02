# -*- coding: utf-8 -*-
"""
Requires numpy, scipy, astropy, matplotlib
@author: A. G. Sreejith
"""
import cute_snr_calculator as csc

t_star = 7000            # Stellar Temperature in K (required)
r_star = 1.2             # Stellar Radius in R sun ()
m_star = 7               # Stellar V magnitude ()          
s_dist = 10              # Stellar parallax in milli arc second
Ra     = 130             # RA in degrees
Dec    = 20              # Dec in degrees   
fwhm=0.8                 # CUTE wavelength resolution in Amstrong
r_noise=3.6              # Readout noise in e-/pix
dark_noise=0.012         # Dark noise in e-/pix/s
exptime=300              # Exposure time in seconds
G=1.0                    # Gain of the CCD
width=10                 # Spectrum width for extraction in pixels
line_core_emission = 1   # Set to one to add MgII lien core emission
mg2_col = None           # MgII column density (will be calculated if set to None)
mg1_col = None           # MgI column density (will be calculated if set to None)
fe2_col = None           # FeII column density (will be calculated if set to None)
add_ism_abs= 1           # Set to one to add ISM absorption
transit_duration=1       # Transit duration in hours   
readtime=20              # CCD read time in seconds
number_of_transits=10    # Number of transits observed
spectral_type='A0V'      # Stellar spectral type
logr=-5.0                # Stellar activity parameter (Log R'HK)
ud=[[2600,2700],[2650,2750]]

fig,snr_txt=csc.cute_snr_calculator(t_star,r_star,m_star,s_dist,logr,spectral_type,Ra,Dec,
                                    transit_duration,ud,fwhm,r_noise,dark_noise,exptime,G,
                                    width,line_core_emission,mg2_col, mg1_col ,fe2_col,
                                    add_ism_abs,readtime,number_of_transits)
filename='cute_snr_'
imagefile=filename+'.png'
txtfile=filename+'.txt'
plt.savefig(os.path.join(os.path.curdir,imagefile),dpi=100)
file_out = os.path.join(os.path.curdir,'static',txtfile)
cute_snr = open(file_out,'w')
cute_snr.write(snr_txt)
cute_snr.close()

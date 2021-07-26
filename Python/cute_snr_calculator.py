"""


@author: A. G. Sreejith
"""


#########################################
###    Import Libraries and Functions
#########################################
import os
import numpy as np
import matplotlib.pyplot as plt
#from PyAstronomy import pyasl
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii
import cute_snr_flux as csf
import extinction as extinction
import unred as unred
import csc_functions as csc
from matplotlib.gridspec import GridSpec


def cute_snr_calculator(t_star,r_star,s_dist,logr,Ra,Dec,transit_duration,td,ud,
                        fwhm=0.8,r_noise=3.6,dark_noise=0.012,exptime=300,
                        gain=1.0,width=10,line_core_emission=1,mg2_col=None,
                        mg1_col=None,fe2_col=None,add_ism_abs=1,readtime=20,
                        number_of_transits=1):
    
    #Constants
    R_sun=6.957e10          #in cm
    
    data=ascii.read(os.path.join(os.path.curdir,'extra','stellar_param_mamjeck.txt')
                    ,comment='#')
        
    Sp=data['col1']
    T_book=data['col2']
    R_book=data['col9']
    Teff=float(t_star)
    loc=csc.find_nearest(T_book,Teff)
    try:
        stype
    except NameError:
        stype  = Sp[loc] 
    try:
        r_star
    except NameError:
        r_star  = R_book[loc]
    try:
        fwhm
    except NameError:
        fwhm = 0.8
    try:
        r_noise
    except NameError:
        r_noise = 3.6
    try:
        dark_noise
    except NameError:
        dark_noise = 0.012
    try:
        exptime
    except NameError:
        exptime = 300
    try:
        gain
    except NameError:
        G = 1
    try:
        width
    except NameError:
        width = 10
    try:
        line_core_emission
    except NameError:    
        line_core_emission = 0
    try:
        mg2_col
    except NameError:
        mg2_col = None
    try:
        mg1_col
    except NameError:
        mg1_col = None
    try:
        fe2_col
    except NameError:
        fe2_col = None
    try:
        add_ism_abs
    except NameError:    
        add_ism_abs = 0
    try:
        readtime
    except NameError:    
        readtime = 20
    try:
        number_of_transits
    except NameError:    
        number_of_transits = 1    
    
    G=float(gain)
    r_star     = r_star*R_sun
    filename   = 't'+str(int(round(t_star,-2))).zfill(5)+'g4.4/model.flx'
    file       = os.path.join(os.path.curdir,'models', filename)
    if os.path.isfile(file):
        fdata  = np.loadtxt(file)
    else:
        filename   = 't'+str(int(round(t_star,-2))+100).zfill(5)+'g4.4/model.flx'
        file       = os.path.join(os.path.curdir,'models', filename)
        fdata      = np.loadtxt(file)
    flux1      = np.zeros(np.shape(fdata))
    flux       = np.zeros(np.shape(fdata))
    flux1[:,1] = (3e18*fdata[:,1])/(fdata[:,0]*fdata[:,0])    #convert to ergs/cm2/s/A
    flux1[:,2] = (3e18*fdata[:,2])/(fdata[:,0]*fdata[:,0])    #convert to ergs/cm2/s/A
    flux[:,0]  = fdata[:,0]
    flux[:,1]  = flux1[:,1]*4*np.pi*(r_star**2)*4*np.pi   #convert to ergs/s/A second 4*!pi for steradian conversion
    flux[:,2]  = flux1[:,2]*4*np.pi*(r_star**2)*4*np.pi   #convert to ergs/s/A second 4*!pi for steradian conversion
    if line_core_emission == 1: 
        data   = csf.cute_snr_lca(flux,0.257,0.288,t_star,r_star,logr,stype) 
    else:
        data   = flux 
    wave       = fdata[:,0]  
    pax        = float(s_dist)                      #paralax in milliarcsec
    d          = 1000.0/pax
    dk         = d/1000.0
    #get extinction based on coordinates and distance
    c      = SkyCoord(ra=Ra, dec=Dec, unit=(u.degree, u.degree))
    glon   = c.galactic.l.deg
    glat   = c.galactic.b.deg
    ebv,av = extinction.extinction_amores(glon,glat,dk) 
    if add_ism_abs == 1:
        if mg2_col == None:
            nh          = 5.8e21*ebv  #The Mg2 column density is
            fractionMg2 = 0.825       #(Frisch & Slavin 2003; this is the fraction of Mg in the ISM that is singly ionised)
            Mg_abn      = -5.33       #(Frisch & Slavin 2003; this is the ISM abundance of Mg)
            nmg2        = np.log10(nh*fractionMg2*10.**Mg_abn)
        if mg1_col == None:
            nh          = 5.8e21*ebv  #The Mg1 column density is
            fractionMg1 = 0.00214     #(Frisch & Slavin 2003; this is the fraction of Mg in the ISM that is singly ionised)
            Mg_abn      = -5.33 ;     #(Frisch & Slavin 2003; this is the ISM abundance of Mg)
            nmg1        = np.log10(nh*fractionMg1*10.**Mg_abn)
        if fe2_col == None:
            nh          = 5.8e21*ebv  #The Fe2 column density is
            fractionFe2 = 0.967       #(Frisch & Slavin 2003; this is the fraction of Fe in the ISM that is singly ionised)
            Fe_abn      = -5.73       #(Frisch & Slavin 2003; this is the ISM abundance of Fe)
            nfe2        = np.log10(nh*fractionFe2*10.**Fe_abn)
         
        flux_data = csf.cute_ism_abs_all(data,nmg2,nmg1,nfe2) 
    else:
        flux_data = data
    flux_di   = flux_data[:,1]
    flux_e    = flux_di/(4.*np.pi*(d*(3.086e+18))**2) #flux at earth
    ebv       = -1.*ebv
    flux_n = unred.unred(wave, flux_e, ebv=ebv, R_V=3.1)
    #flux_n=flux_e
    #Useful defs
    #    1 Jy = 10^-23 erg sec^-1 cm^-2 Hz^-1
    #    1 Jy = 1.51e7 photons sec^-1 m^-2 (dlambda/lambda)^-1
    #convert to photons
    photons_star    = flux_n*5.03e7*wave     #from ergs/s/cm2/A to photons/s/cm2/A
    #photons_star    = np.zeros(2,len(wave))
    st              = csc.find_nearest(wave, 2000)
    en              = csc.find_nearest(wave, 4000)
    wave_new        = wave[st:en]
    photons_star_new= photons_star[st:en]
    #convolution with instrument response
    smoothedflux    = csc.gaussbroad(wave_new,photons_star_new,fwhm/2.0)
    file_wave       = os.path.join(os.path.curdir,'extra', 'wavelength.txt')
    eff_file        = os.path.join(os.path.curdir,'extra', 'eff_area.txt')
    wavelength      = np.loadtxt(file_wave)
    w_length        = len(wavelength)
    ccd_flux        = np.zeros(w_length)
    ccd_wave        = np.zeros(w_length)
    wave_res        = np.zeros(int(w_length/2))
    for i in range(0,w_length-1,2):
        j=int(i/2)
        wave_res[j] = (wavelength[i]+wavelength[i+1])/2
    #interpolate and trim to detector size
    ccd_flux  = np.interp(wave_res,wave_new,smoothedflux)
    eff_area  = np.loadtxt(eff_file)
    aeff      = np.interp(wave_res,eff_area[:,0],eff_area[:,4])
    ccd_count1= np.zeros(int(w_length/2))
    ccd_count1= ccd_flux*aeff*fwhm 
    ccd_count = np.zeros(w_length)
    noise     = np.zeros(w_length)
    snr       = np.zeros(w_length)   
    
    for i in range (0,w_length):
        ccd_count[i]  = (ccd_count1[int(i/2)]/2)*exptime*G # assuming 2 resolution element
        noise[i]      = np.sqrt(ccd_count[i]+(width*(r_noise**2+(dark_noise*exptime*G))))
        snr[i]        = ccd_count[i]/noise[i]

#user defined region
    user_low=[item[0] for item in ud]
    user_high=[item[1] for item in ud]
    #SNR calculations
    nx            = w_length
    spectra_1dwf  = ccd_count
    error_1dwf    = noise
    lenby3        = int(nx/3)
    st            = 1
    en            = nx
    n_w           = en-st+1
    x1            = wavelength[st:en]
    y1            = spectra_1dwf[st:en]
    dy1           = error_1dwf[st:en]
    tpf,tpe       = csc.trapz_error(x1,y1,dy1)
    snr_full      = tpf/tpe
    st            = 0
    en            = lenby3
    n_w           = en-st
    x2            = wavelength[st:en]
    y2            = spectra_1dwf[st:en]
    dy2           = error_1dwf[st:en]
    tpf,tpe       = csc.trapz_error(x2,y2,dy2)
    snr_blue      = tpf/tpe
    st            = lenby3
    en            = 2*lenby3
    n_w           = en-st
    x3            = wavelength[st:en]
    y3            = spectra_1dwf[st:en]
    dy3           = error_1dwf[st:en]
    tpf,tpe       = csc.trapz_error(x3,y3,dy3)
    snr_middle    = tpf/tpe
    st            = 2*lenby3
    en            = nx
    x4            = wavelength[st:en]
    y4            = spectra_1dwf[st:en]
    dy4           = error_1dwf[st:en]
    tp,te         = csc.trapz_error(x4,y4,dy4)
    snr_red       = tpf/tpe
    
    # MgII 2795,2802
    st            = csc.find_nearest(wavelength, 2793)
    en            = csc.find_nearest(wavelength, 2805)
    x5            = wavelength[st:en]
    y5            = spectra_1dwf[st:en]
    dy5           = error_1dwf[st:en]
    tpf,tpe       = csc.trapz_error(x5,y5,dy5)
   
    snr_mg2       = tpf/tpe
    #MgI 2852
    st            = csc.find_nearest(wavelength, 2850)
    en            = csc.find_nearest(wavelength, 2854)
    x6            = wavelength[st:en]
    y6            = spectra_1dwf[st:en]
    dy6           = error_1dwf[st:en]
    tpf,tpe       = csc.trapz_error(x6,y6,dy6)
    snr_mg1       = tpf/tpe
    #Fe II 2585
    st            = csc.find_nearest(wavelength, 2583)
    en            = csc.find_nearest(wavelength, 2587)
    x7            = wavelength[st:en]
    y7            = spectra_1dwf[st:en]
    dy7           = error_1dwf[st:en]
    tpf,tpe       = csc.trapz_error(x7,y7,dy7)
    snr_fe2       = tpf/tpe
    
    #number of exposures in transit duration
    t_dur         = transit_duration*3600. #in seconds
    obs_time      = exptime+int(readtime)  #in seconds
    obs_transit   = t_dur/obs_time         #number of observations in transit
    
    #transit depth calculations per transit
    unc_transit_full   = np.sqrt(2.)/(np.sqrt(obs_transit)*snr_full)
    unc_transit_blue   = np.sqrt(2.)/(np.sqrt(obs_transit)*snr_blue)
    unc_transit_middle = np.sqrt(2.)/(np.sqrt(obs_transit)*snr_middle)
    unc_transit_red    = np.sqrt(2.)/(np.sqrt(obs_transit)*snr_red)
    unc_transit_mg1    = np.sqrt(2.)/(np.sqrt(obs_transit)*snr_mg1)
    unc_transit_mg2    = np.sqrt(2.)/(np.sqrt(obs_transit)*snr_mg2)
    unc_transit_fe2    = np.sqrt(2.)/(np.sqrt(obs_transit)*snr_fe2)
    
    depth = np.sqrt((1.0-td)/100)
    unc_radius_full    = np.sqrt(depth)/(np.sqrt(2.)*np.sqrt(obs_transit)*depth*snr_full)
    unc_radius_blue    = np.sqrt(depth)/(np.sqrt(2.)*np.sqrt(obs_transit)*depth*snr_blue)
    unc_radius_middle  = np.sqrt(depth)/(np.sqrt(2.)*np.sqrt(obs_transit)*depth*snr_middle)
    unc_radius_red     = np.sqrt(depth)/(np.sqrt(2.)*np.sqrt(obs_transit)*depth*snr_red)
    unc_radius_mg1     = np.sqrt(depth)/(np.sqrt(2.)*np.sqrt(obs_transit)*depth*snr_mg1)
    unc_radius_mg2     = np.sqrt(depth)/(np.sqrt(2.)*np.sqrt(obs_transit)*depth*snr_mg2)
    unc_radius_fe2     = np.sqrt(depth)/(np.sqrt(2.)*np.sqrt(obs_transit)*depth*snr_fe2)
    
    #transit depth calculations for n transit
    n=number_of_transits
    n_unc_transit_full   = unc_transit_full/np.sqrt(n)
    n_unc_transit_blue   = unc_transit_blue/np.sqrt(n)
    n_unc_transit_middle = unc_transit_middle/np.sqrt(n)
    n_unc_transit_red    = unc_transit_red/np.sqrt(n)
    n_unc_transit_mg1    = unc_transit_mg1/np.sqrt(n)
    n_unc_transit_mg2    = unc_transit_mg2/np.sqrt(n)
    n_unc_transit_fe2    = unc_transit_fe2/np.sqrt(n)
    
    cute_snr='-------------------------------------------------------------------------------\n'
    cute_snr=cute_snr+'                      CUTE SNR CALCULATIONS\n'
    cute_snr=cute_snr+'-------------------------------------------------------------------------------\n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Input Parameters \n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Stellar Parameters \n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Stellar Temperature(K)  : '+str(t_star)+'\n'
    cute_snr=cute_snr+'Stellar Radius (R_sun)  : '+str(r_star/R_sun)+'\n'
    #cute_snr=cute_snr+'Stellar Magnitude (V)   : '+str(m_star)+'\n'
    cute_snr=cute_snr+'Stellar Distance (mas)  : '+str(s_dist)+'\n'
    cute_snr=cute_snr+'Target RA (deg)         : '+str(Ra)+'\n'
    cute_snr=cute_snr+'Target Declination (deg): '+str(Dec)+'\n'
    cute_snr=cute_snr+'\n'
    if line_core_emission == 1:
        cute_snr=cute_snr+'Line core emission added.\n'
        cute_snr=cute_snr+'Stellar Specrtral type  : '+str(stype)+'\n'
        cute_snr=cute_snr+"logR'HK                 : "+str(logr)+'\n' 
    if add_ism_abs == 1:
        cute_snr=cute_snr+'ISM absorption added for the following species \n'
        if mg1_col == 1: 
            cute_snr=cute_snr+'MgI column density      : '+str(mg1_col)+'\n'
        else:
            cute_snr=cute_snr+'MgI column density was calculated.\n'
        if mg2_col == 1:
            cute_snr=cute_snr+'MgII column density     : '+str(mg2_col)+'\n'
        else:
            cute_snr=cute_snr+'MgII column density was calculated. \n'
        if fe2_col == 1:  
            cute_snr=cute_snr+'FeII column density     : '+str(fe2_col)+'\n'
        else:
            cute_snr=cute_snr+'FeII column density was calculated.\n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Instrument Parameters \n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Spectral Resolution (A) : '+str(fwhm)+'\n'
    cute_snr=cute_snr+'Spectrum Width (pix)    : '+str(width)+'\n'
    cute_snr=cute_snr+'in cross-disp direction \n'
    cute_snr=cute_snr+'Readout Noise (e/pix)   : '+str(r_noise)+'\n'
    cute_snr=cute_snr+'Dark Noise (e/pix/s)    : '+str(dark_noise)+'\n' 
    cute_snr=cute_snr+'Exposure Time (s)       : '+str(exptime)+'\n'
    cute_snr=cute_snr+'Read Time (s)           : '+str(readtime)+'\n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Transit Parameters \n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Transit Duration (hrs)  : '+str(transit_duration)+'\n'
    cute_snr=cute_snr+'Number of Transits      : '+str(number_of_transits)+'\n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Input Files \n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Stellar Model file      : LL model:'+file+'\n'
    cute_snr=cute_snr+'Wavelength file         : '+str(file_wave)+'\n'
    cute_snr=cute_snr+'Effective area file     : '+str(eff_file)+'\n'
    cute_snr=cute_snr+'-------------------------------------------------------------------------------\n'
    cute_snr=cute_snr+'                        CUTE SNR OUTPUTS\n'
    cute_snr=cute_snr+'-------------------------------------------------------------------------------\n'
    cute_snr=cute_snr+'Wavelength Region [A]\t\t\t SNR\t Uncertainty in Transit Depth [ppm]\n'
    cute_snr=cute_snr+'                                             1 Transit\t'+str(n)+' Transits\n'  
    cute_snr=cute_snr+'-------------------------------------------------------------------------------\n'  
    cute_snr=cute_snr+'Full Band ['+str(round(x1[0],2))+'-'+str(round(x1[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(snr_full,4))+'\t'+str(round(unc_transit_full*1E6,4))+'\t'
    cute_snr=cute_snr+str(round(n_unc_transit_full*1E6,4))+'\n'
    cute_snr=cute_snr+'Lower Band ['+str(round(x2[0],2))+'-'+str(round(x2[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(snr_blue,4))+'\t'+str(round(unc_transit_blue*1E6,4))+'\t'
    cute_snr=cute_snr+str(round(n_unc_transit_blue*1E6,4))+'\n'
    cute_snr=cute_snr+'Mid Band ['+str(round(x3[0],2))+'-'+str(round(x3[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(snr_middle,4))+'\t'+str(round(unc_transit_middle*1E6,4))+'\t'
    cute_snr=cute_snr+str(round(n_unc_transit_middle*1E6,4))+'\n'
    cute_snr=cute_snr+'Upper Band ['+str(round(x4[0],2))+'-'+str(round(x4[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(snr_red,4))+'\t'+str(round(unc_transit_red*1E6,4))+'\t'
    cute_snr=cute_snr+str(round(n_unc_transit_red*1E6,4))+'\n'
    cute_snr=cute_snr+'MgII Band ['+str(round(x5[0],2))+'-'+str(round(x5[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(snr_mg2,4))+'\t'+str(round(unc_transit_mg2*1E6,4))+'\t'
    cute_snr=cute_snr+str(round(n_unc_transit_mg2*1E6,4))+'\n'
    cute_snr=cute_snr+'MgI Band ['+str(round(x6[0],2))+'-'+str(round(x6[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(snr_mg1,4))+'\t'+str(round(unc_transit_mg1*1E6,4))+'\t'
    cute_snr=cute_snr+str(round(n_unc_transit_mg1*1E6,4))+'\n'
    cute_snr=cute_snr+'FeII Band ['+str(round(x7[0],2))+'-'+str(round(x7[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(snr_fe2,4))+'\t'+str(round(unc_transit_fe2*1E6,4))+'\t'
    cute_snr=cute_snr+str(round(n_unc_transit_fe2*1E6,4))+'\n'
                   
    for i in range(0,len(user_low)):
        ud_low    = user_low[i]
        ud_high   = user_high[i]
        st        = csc.find_nearest(wavelength, ud_low)
        en        = csc.find_nearest(wavelength, ud_high)
        x_ud      = wavelength[st:en]
        y_ud      = spectra_1dwf[st:en]
        dy_ud     = error_1dwf[st:en]
        tpf,tpe   = csc.trapz_error(x_ud,y_ud,dy_ud)
        snr_ud = tpf/tpe
        unc_transit_ud    = 1./(np.sqrt(obs_transit)*snr_ud)
        n_unc_transit_ud   = unc_transit_ud/np.sqrt(n)
        cute_snr=cute_snr+'User Band '+str(i+1)+' ['+str(round(x_ud[0],2))+'-'+str(round(x_ud[-1],2))+']\t\t'
        cute_snr=cute_snr+str(round(snr_ud,4))+'\t'+str(round(unc_transit_ud*1E6,4))+'\t'
        cute_snr=cute_snr+str(round(n_unc_transit_ud*1E6,4))+'\n' 
    cute_snr=cute_snr+'-------------------------------------------------------------------------------\n'
    cute_snr=cute_snr+'Radius uncertanitiy\n'
    cute_snr=cute_snr+'Wavelength Region [A]\t\t\t Uncertainty in Radius [ppm] for 1 Transits\n'  
    cute_snr=cute_snr+'-------------------------------------------------------------------------------\n'  
    cute_snr=cute_snr+'Full Band ['+str(round(x1[0],2))+'-'+str(round(x1[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(unc_radius_full*1E6,4))+'\n'
    cute_snr=cute_snr+'Lower Band ['+str(round(x2[0],2))+'-'+str(round(x2[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(unc_radius_blue*1E6,4))+'\n'
    cute_snr=cute_snr+'Mid Band ['+str(round(x3[0],2))+'-'+str(round(x3[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(unc_radius_middle*1E6,4))+'\n'
    cute_snr=cute_snr+'Upper Band ['+str(round(x4[0],2))+'-'+str(round(x4[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(unc_radius_red*1E6,4))+'\n'
    cute_snr=cute_snr+'MgII Band ['+str(round(x5[0],2))+'-'+str(round(x5[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(unc_radius_mg2*1E6,4))+'\n'
    cute_snr=cute_snr+'MgI Band ['+str(round(x6[0],2))+'-'+str(round(x6[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(unc_radius_mg1*1E6,4))+'\n'
    cute_snr=cute_snr+'FeII Band ['+str(round(x7[0],2))+'-'+str(round(x7[-1],2))+']\t\t'
    cute_snr=cute_snr+str(round(unc_radius_fe2*1E6,4))+'\n'
                   
    for i in range(0,len(user_low)):
        ud_low    = user_low[i]
        ud_high   = user_high[i]
        st        = csc.find_nearest(wavelength, ud_low)
        en        = csc.find_nearest(wavelength, ud_high)
        x_ud      = wavelength[st:en]
        y_ud      = spectra_1dwf[st:en]
        dy_ud     = error_1dwf[st:en]
        tpf,tpe   = csc.trapz_error(x_ud,y_ud,dy_ud)
        snr_ud = tpf/tpe
        unc_radius_ud    = np.sqrt(depth)/(np.sqrt(2.)*depth*np.sqrt(obs_transit)*snr_ud)
        cute_snr=cute_snr+'User Band '+str(i+1)+' ['+str(round(x_ud[0],2))+'-'+str(round(x_ud[-1],2))+']\t\t'
        cute_snr=cute_snr+str(round(unc_radius_ud*1E6,4))+'\n' 
    cute_snr=cute_snr+'-------------------------------------------------------------------------------\n'
    
    cute_snr=cute_snr+'Full Spectrum SNR \n'
    cute_snr=cute_snr+'\n'
    cute_snr=cute_snr+'Wavelength [A]\tFlux [counts] \tNoise[counts] \tSNR \n'
    for i in range(w_length):
        cute_snr=cute_snr+str(round(wavelength[i],2))+'\t\t\t'+str(round(ccd_count[i],4))+'\t\t'
        cute_snr=cute_snr+str(round(noise[i],4))+'\t\t'+str(round(snr[i],4))+'\n'

    fig = plt.figure(constrained_layout=True)
    gs = GridSpec(3, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(wavelength, ccd_count, '-', color='red')
    ax1.fill_between(wavelength, ccd_count - noise, ccd_count + noise,
                 color='lightsteelblue', alpha=0.2)
    ax1.set_xlabel('wavelength[$\AA$]')
    ax1.set_ylabel('flux[counts]')
    ax1.set_title('CUTE full band')
    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(wavelength, snr, '-', color='red')
    ax2.set_xlabel('wavelength[$\AA$]')
    ax2.set_ylabel('SNR')
    ax2.set_title('CUTE full band SNR')
    ax3= fig.add_subplot(gs[2, 0])
    ax3.plot(wavelength[670:740], ccd_count[670:740], '-', color='red')
    ax3.fill_between(wavelength[670:740], ccd_count[670:740] - noise[670:740],
                     ccd_count[670:740] + noise[670:740],
                     color='lightsteelblue', alpha=0.2)
    ax3.set_xlabel('wavelength[$\AA$]')
    ax3.set_ylabel('flux[counts]')
    ax3.set_title('CUTE MgII')
    ax4= fig.add_subplot(gs[2, 1])
    ax4.plot(wavelength[1213:1464], ccd_count[1213:1464], '-', color='red')
    ax4.fill_between(wavelength[1213:1464], ccd_count[1213:1464] - noise[1213:1464], 
                     ccd_count[1213:1464] + noise[1213:1464],
                 color='lightsteelblue', alpha=0.2)
    ax4.set_xlabel('wavelength[$\AA$]')
    ax4.set_ylabel('flux[counts]')
    ax4.set_title('CUTE Continuum')
    fig.show()
    #    plt.savefig(os.path.join(os.path.curdir,'output','cute.png'),dpi=100)
    return fig,cute_snr
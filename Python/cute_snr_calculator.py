# -*- coding: utf-8 -*-
"""


@author: A. G. Sreejith
"""


#########################################
###    Import Libraries and Functions
#########################################
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii
import cute_snr_flux as csf
import extinction as extinction
import csc_functions as csc


def cute_snr_calculator(t_star,r_star,m_star,s_dist,logr,stype,Ra,Dec,transit_duration,ud,fwhm=0.8,r_noise=3.6,dark_noise=0.012,
                        exptime=300,G=1,width=10,line_core_emission=1,mg2_col=None,mg1_col=None,fe2_col=None,
                        add_ism_abs=1,readtime=20,number_of_transits=1):
    
    print('SNR begins')
    #Constants
    R_sun=6.957e10          #in cm
    L_sun=3.828e33
    
    data=ascii.read(os.path.join(os.path.curdir,'extra','stellar_param_mamjeck.txt'),comment='#')
        
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
    t          = float(t_star)
    t4         = t**4
    r2         = r_star**2
    stepahs    = 5.6704e-5
    L          = stepahs*4*np.pi*r2*t4
    if line_core_emission == 1: 
        data   = csf.cute_snr_lca(flux,0.257,0.288,t_star,r_star,logr,stype) 
    else:
        data   = flux 
    wave       = fdata[:,0]  
    flux       = fdata[:,1]  
    cont       = fdata[:,2]
    pax        = float(s_dist)                      #paralax in milliarcsec
    d          = 1000.0/pax
    dk         = d/1000.0
    #get extinction based on coordinates and distance
    c      = SkyCoord(ra=Ra, dec=Dec, unit=(u.degree, u.degree))
    glon   = c.galactic.l.deg
    glat   = c.galactic.b.deg
    ebv,av = extinction.extinction_amores(glon,glat,dk) 
    
    
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
    if add_ism_abs == 1:
        flux_data = csf.cute_ism_abs_all(data,nmg2,nmg1,nfe2) 
    else:
        flux_data = data
    flux_di   = flux_data[:,1]
    flux_e    = flux_di/(4.*np.pi*(d*(3.086e+18))**2) #flux at earth
    ebv       = -1.*ebv
    #flux_n = pyasl.unred(wave, flux_e, ebv=ebv, R_V=3.1)
    flux_n=flux_e
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
        ccd_count[i]  = (ccd_count1[int(i/2)]/2)*exptime # assuming 2 resolution element 
        noise[i]      = np.sqrt(ccd_count[i]+(width*(r_noise**2+(dark_noise*exptime))))
        snr[i]        = ccd_count[i]/noise[i]

#user defined region
    user_low=ud[0][:]
    user_high=ud[1][:]

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
    n_w           = en-st+1
    x4            = wavelength[st:en]
    y4            = spectra_1dwf[st:en]
    dy4           = error_1dwf[st:en]
    tp,te         = csc.trapz_error(x4,y4,dy4)
    snr_red       = tpf/tpe
    
    # MgII 2795,2802
    st            = csc.find_nearest(wavelength, 2793)
    en            = csc.find_nearest(wavelength, 2805)
    n_w           = en-st
    x5            = wavelength[st:en]
    y5            = spectra_1dwf[st:en]
    dy5           = error_1dwf[st:en]
    tpf,tpe       = csc.trapz_error(x5,y5,dy5)
   
    snr_mg2       = tpf/tpe
    #MgI 2852
    st            = csc.find_nearest(wavelength, 2850)
    en            = csc.find_nearest(wavelength, 2854)
    n_w           = en-st+1
    x6            = wavelength[st:en]
    y6            = spectra_1dwf[st:en]
    dy6           = error_1dwf[st:en]
    tpf,tpe       = csc.trapz_error(x6,y6,dy6)
    snr_mg1       = tpf/tpe
    #Fe II 2585
    st            = csc.find_nearest(wavelength, 2583)
    en            = csc.find_nearest(wavelength, 2587)
    n_w           = en-st
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
    unc_tranit_full   = 1./(np.sqrt(obs_transit)*snr_full)
    unc_tranit_blue   = 1./(np.sqrt(obs_transit)*snr_blue)
    unc_tranit_middle = 1./(np.sqrt(obs_transit)*snr_middle)
    unc_tranit_red    = 1./(np.sqrt(obs_transit)*snr_red)
    unc_tranit_mg1    = 1./(np.sqrt(obs_transit)*snr_mg1)
    unc_tranit_mg2    = 1./(np.sqrt(obs_transit)*snr_mg2)
    unc_tranit_fe2    = 1./(np.sqrt(obs_transit)*snr_fe2)
    
    #transit depth calculations for n transit
    n=number_of_transits
    n_unc_tranit_full   = unc_tranit_full/np.sqrt(n)
    n_unc_tranit_blue   = unc_tranit_blue/np.sqrt(n)
    n_unc_tranit_middle = unc_tranit_middle/np.sqrt(n)
    n_unc_tranit_red    = unc_tranit_red/np.sqrt(n)
    n_unc_tranit_mg1    = unc_tranit_mg1/np.sqrt(n)
    n_unc_tranit_mg2    = unc_tranit_mg2/np.sqrt(n)
    n_unc_tranit_fe2    = unc_tranit_fe2/np.sqrt(n)
    
    file_out = os.path.join(os.path.curdir,'output','cute_snr.txt')
    cute_snr = open(file_out,'w')
    cute_snr.write('-------------------------------------------------------------------------------\n'
                   +'                      CUTE SNR CALCULATIONS\n'
                   +'-------------------------------------------------------------------------------\n'
                   +'\n'
                   +'Input Parameters \n'
                   +'\n'
                   +'Stellar Parameters \n'
                   +'\n'
                   +'Stellar Temperature(K)  : '+str(t_star)+'\n'
                   +'Stellar Radius (R_sun)  : '+str(r_star/R_sun)+'\n'
                   +'Stellar Magnitude (V)   : '+str(m_star)+'\n'
                   +'Stellar Distance (mas)  : '+str(s_dist)+'\n'
                   +'Stellar Specrtral type  : '+str(stype)+'\n'
                   +'Target RA (deg)         : '+str(Ra)+'\n'
                   +'Target Declination (deg): '+str(Dec)+'\n'
                   +'\n')
    if line_core_emission == 1:
        cute_snr.write('Line core emission added.\n'
                       +"logR'HK                 : "+str(logr)+'\n') 
    if add_ism_abs == 1:
        cute_snr.write('ISM absorption added for the following species \n')
        if mg1_col == 1: 
            cute_snr.write('MgI column density      : '+str(mg1_col)+'\n')
        else:
            cute_snr.write('MgI column density was calculated.\n')
        if mg2_col == 1:
            cute_snr.write('MgII column density     : '+str(mg2_col)+'\n')
        else:
            cute_snr.write('MgII column density was calculated. \n')
        if fe2_col == 1:  
            cute_snr.write('FeII column density     : '+str(fe2_col)+'\n')
        else:
            cute_snr.write('FeII column density was calculated.\n')
    cute_snr.write('\n'
                   +'Instrument Parameters \n'
                   +'\n'
                   +'Spectral Resolution (A) : '+str(fwhm)+'\n'
                   +'Spectrum Width (pix)    : '+str(width)+'\n'
                   +'in cross-disp direction \n'
                   +'Readout Noise (e/pix)   : '+str(r_noise)+'\n'
                   +'Dark Noise (e/pix/s)    : '+str(dark_noise)+'\n' 
                   +'Exposure Time (s)       : '+str(exptime)+'\n'
                   +'Read Time (s)           : '+str(readtime)+'\n'
                   +'\n'
                   +'Input Files \n'
                   +'\n'
                   +'Stellar Model file      : LL model:'+file+'\n'
                   +'Wavelength file         : '+str(file_wave)+'\n'
                   +'Effective area file     : '+str(eff_file)+'\n'
                   +'-------------------------------------------------------------------------------\n'
                   +'                        CUTE SNR OUTPUTS\n'
                   +'-------------------------------------------------------------------------------\n'
                   +'Wavelength Region [A]\t\t\t SNR\t Uncertanity in Transit Depth [ppm]\n'
                   +'                                             1 Transit\t'+str(n)+' Transits\n'  
                   +'-------------------------------------------------------------------------------\n'  
                   +'Full Band['+str(round(x1[0],2))+'-'+str(round(x1[-1],2))+']\t\t'+
                   str(round(snr_full,4))+'\t'+str(round(unc_tranit_full*1E6,4))+'\t'
                   +str(round(n_unc_tranit_full*1E6,4))+'\n'
                   +'Lower Band['+str(round(x2[0],2))+'-'+str(round(x2[-1],2))+']\t\t'
                   +str(round(snr_blue,4))+'\t'+str(round(unc_tranit_blue*1E6,4))+'\t'
                   +str(round(n_unc_tranit_blue*1E6,4))+'\n'
                   +'Mid Band['+str(round(x3[0],2))+'-'+str(round(x3[-1],2))+']\t\t'
                   +str(round(snr_middle,4))+'\t'+str(round(unc_tranit_middle*1E6,4))+'\t'
                   +str(round(n_unc_tranit_middle*1E6,4))+'\n'
                   +'Upper Band['+str(round(x4[0],2))+'-'+str(round(x4[-1],2))+']\t\t'
                   +str(round(snr_red,4))+'\t'+str(round(unc_tranit_red*1E6,4))+'\t'
                   +str(round(n_unc_tranit_red*1E6,4))+'\n'
                   +'MgII Band['+str(round(x5[0],2))+'-'+str(round(x5[-1],2))+']\t\t'+
                   str(round(snr_mg2,4))+'\t'+str(round(unc_tranit_mg2*1E6,4))+'\t'+
                   str(round(n_unc_tranit_mg2*1E6,4))+'\n'
                   +'MgI Band['+str(round(x6[0],2))+'-'+str(round(x6[-1],2))+']\t\t'+
                   str(round(snr_mg1,4))+'\t'+str(round(unc_tranit_mg1*1E6,4))+'\t'+
                   str(round(n_unc_tranit_mg1*1E6,4))+'\n'
                   +'FeII Band['+str(round(x7[0],2))+'-'+str(round(x7[-1],2))+']\t\t'+
                   str(round(snr_fe2,4))+'\t'+str(round(unc_tranit_fe2*1E6,4))+'\t'+
                   str(round(n_unc_tranit_fe2*1E6,4))+'\n') 
    for i in range(0,len(user_low)):
        ud_low    = user_low[i]
        ud_high   = user_high[i]
        st        = csc.find_nearest(wavelength, ud_low)
        en        = csc.find_nearest(wavelength, ud_high)
        n_w       = en-st+1
        x_ud      = wavelength[st:en]
        y_ud      = spectra_1dwf[st:en]
        dy_ud     = error_1dwf[st:en]
        tpf,tpe   = csc.trapz_error(x_ud,y_ud,dy_ud)
        snr_ud = tpf/tpe
        unc_tranit_ud    = 1./(np.sqrt(obs_transit)*snr_ud)
        n_unc_tranit_ud   = unc_tranit_ud/np.sqrt(n)
        cute_snr.write('User Band '+str(i+1)+'['+str(round(x_ud[0],2))+'-'+str(round(x_ud[-1],2))+']\t\t'+
                       str(round(snr_ud,4))+'\t'+str(round(unc_tranit_ud*1E6,4))+'\t'+
                       str(round(n_unc_tranit_ud*1E6,4))+'\n')  
    cute_snr.write('-------------------------------------------------------------------------------\n'
               +'Full Spectrum SNR \n'
               +'\n'
               +'Wavelength [A]\tFlux [counts] \tNoise[counts] \tSNR \n')
    for i in range(w_length):
        cute_snr.write(str(round(wavelength[i],2))+'\t\t\t'+str(round(ccd_count[i],4))+'\t\t'+
                       str(round(noise[i],4))+'\t\t'+str(round(snr[i],4))+'\n')
    cute_snr.close()

    plt.figure()
    
    plt.plot(wavelength, ccd_count, '-', color='black')
    plt.fill_between(wavelength, ccd_count - noise, ccd_count + noise,
                 color='gray', alpha=0.2)
    plt.xlabel('wavelength[$\AA$]')
    plt.ylabel('counts')
    plt.savefig(os.path.join(os.path.curdir,'output','cute.png'))
 
''''

  cgDisplay, 1600, 1000

  ;POSITION=[leftAxisLoc, bottomAxisLoc, rightAxesLoc, topAxesLoc]
  ;new_wave= make_array((wavelength[w_length-1]-wavelength[0])*10,START=wavelength[0],INCREMENT=0.1,/DOUBLE, /INDEX)
  ;new_count=interpol(ccd_count,wavelength,new_wave,/SPLINE)
  ;new_noise=interpol(noise,wavelength,new_wave,/SPLINE)
  high_error = (ccd_count + noise)
  low_error = (ccd_count - noise)
  ; Draw the line plot with no data
  cgPlot, wavelength, ccd_count, Title='CUTE FULL BAND', XTitle='wavelength [$\Angstrom$]',/NoErase $
      , YTitle='flux [counts]', Position=[0.1, 0.2875, 0.925, 0.475], YStyle=1, /NoData $
      , xrange=[wavelength[0]-10,wavelength[w_length-1]+10], yrange=[min(low_error)-100,max(high_error)+100]

  ; Fill in the error estimates.
  cgColorFill, [wavelength, Reverse(wavelength), wavelength[0]], $
    [high_error, Reverse(low_error), high_error[0]], $
    Color='sky blue'

  ; Draw the line plot with no data
  cgPlotS, wavelength, ccd_count, Color='red';,PSym=-16, SymColor='olive', $
  ;SymSize=1.0, Thick=2
  cgPlot, wavelength, snr, Color='red',XTitle='wavelength [$\Angstrom$]',$
    YTitle='SNR', Position=[0.1, 0.05, 0.925, 0.2375],/NoErase $
      ,xrange=[wavelength[0]-10,wavelength[w_length-1]+10], yrange=[min(snr)-5,max(snr)+5] 


  ;Calculate for region of interest as well
  cgPlot, wavelength[670:740], ccd_count[670:740], Title='CUTE MgII', XTitle='wavelength [$\Angstrom$]',/NoErase $
    , YTitle='flux [counts]', Position=[0.1, 0.7875, 0.475, 0.975], YStyle=1, /NoData $
    , xrange=[wavelength[670]-2,wavelength[740]+2], yrange=[min(low_error[670:740])-50,max(high_error[670:740])+50]

  ; Fill in the error estimates.
  cgColorFill, [wavelength[670:740], Reverse(wavelength[670:740]), wavelength[670]], $
    [high_error[670:740], Reverse(low_error[670:740]), high_error[670]], $
    Color='sky blue'

  ; Draw the line plot with no data
  cgPlotS, wavelength[670:740], ccd_count[670:740], Color='red';,PSym=-16, SymColor='olive', $
  ;SymSize=1.0, Thick=2
  cgPlotS, [2793,2793], !Y.CRange,linestyle=2,color='blue'
  cgPlotS, [2805,2805], !Y.CRange,linestyle=2,color='blue'
  cgPlot, wavelength[670:740], snr[670:740], Color='red',XTitle='wavelength [$\Angstrom$]',$
    YTitle='SNR', Position=[0.1, 0.55, 0.475, 0.7375],/NoErase $
    ,xrange=[wavelength[670]-2,wavelength[740]+2], yrange=[min(snr[670:740])-2,max(snr[670:740])+2]
  cgPlotS, [2793,2793], !Y.CRange,linestyle=2,color='blue'  
  cgPlotS, [2805,2805], !Y.CRange,linestyle=2,color='blue'
  tp=trapz_error(wavelength[692:721],ccd_count[692:721],noise[692:721])
  mgsnr=tp[0]/tp[1]
  cgText,wavelength[720],min(snr[670:740]),'INTEGRATED SNR: '+STRTRIM(string(mgsnr),2),CHARSIZE=1., charthick = 1.,/DATA,color='black'

  ;Continuum region
  cgPlot, wavelength[1213:1464], ccd_count[1213:1464], Title='CUTE Continuum', XTitle='wavelength [$\Angstrom$]',/NoErase $
    , YTitle='flux [counts]', Position=[0.55, 0.7875, 0.925, 0.975], YStyle=1, /NoData $
    , xrange=[wavelength[1213]-2,wavelength[1464]+2], yrange=[min(low_error[1213:1464])-50,max(high_error[1213:1464])+50]

  ; Fill in the error estimates.
  cgColorFill, [wavelength[1213:1464], Reverse(wavelength[1213:1464]), wavelength[1213]], $
    [high_error[1213:1464], Reverse(low_error[1213:1464]), high_error[1213]], $
    Color='sky blue'

  ; Draw the line plot with no data
  cgPlotS, wavelength[1213:1464], ccd_count[1213:1464], Color='red';,PSym=-16, SymColor='olive', $
    ;SymSize=1.0, Thick=2
  ;cgPlotS, [2793,2793], !Y.CRange,linestyle=2,color='blue'
  ;cgPlotS, [2805,2805], !Y.CRange,linestyle=2,color='blue'
  cgPlot, wavelength[1213:1464], snr[1213:1464], Color='red',XTitle='wavelength [$\Angstrom$]',$
    YTitle='SNR', Position=[0.55, 0.55, 0.925, 0.7375],/NoErase $
    ,xrange=[wavelength[1213]-2,wavelength[1464]+2], yrange=[min(snr[1213:1464])-2,max(snr[1213:1464])+2]
  ;cgPlotS, [2793,2793], !Y.CRange,linestyle=2,color='blue'
  ;cgPlotS, [2805,2805], !Y.CRange,linestyle=2,color='blue'
  tp=trapz_error(wavelength[1213:1464],ccd_count[1213:1464],noise[1213:1464])
  mgsnr=tp[0]/tp[1]
  cgText,wavelength[1350],min(snr[1213:1464]-1),'INTEGRATED SNR: '+STRTRIM(string(mgsnr),2),CHARSIZE=1., charthick = 1.,/DATA,color='black'
  write_png,infile.file_out+'cute_snr_plot.png',TVRD(/TRUE)
'''''

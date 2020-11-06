pro cute_snr_calculator,parameter_file

  infile     = gm_read_textstructure(parameter_file)
  file       = infile.input_file
  t_star     = double(infile.stellar_temperature)
  r_star     = double(infile.stellar_radius)
  M_star     = double(infile.stellar_magnitude)
  fwhm       = float(infile.spectral_resolution)
  file_wave  = infile.wave_file
  file_eff   = infile.file_eff
  r_noise    = double(infile.readout_noise)
  exptime    = long(infile.exptime)
  dark_noise = double(infile.dark_noise)
  G          = double(infile.ccd_gain)
  ra         = double(infile.Ra)
  dec        = double(infile.Dec)
  width      = long(infile.spectrum_width)
  ;read input wavelength file
  w_length   = file_lines(file_wave)
  wavelength = dblarr(w_length)
  ccd_flux   = dblarr(w_length)
  ccd_wave   = dblarr(w_length)
  openr,1,file_wave
  readf,1,wavelength
  close,1

  photons           = cute_snr_photons(infile)
  photons_star      = photons[1,*]
  wave              = photons[0,*]
  photons_star      = reform(photons_star,n_elements(photons_star))
  wave              = reform(wave,n_elements(wave))
  st                = value_locate(wave, 2000)
  en                = value_locate(wave, 7000)
  wave_new          = wave[st:en]
  photons_star_new  = photons_star[st:en]
  hwhm              = fwhm/2
  ;convolution with instrument response
  smoothedflux      = gaussbroad(wave_new,photons_star_new,hwhm)
  wave_res          = dblarr(w_length/2)
  for i=0,w_length-2,2 do wave_res[i/2]=mean(wavelength[i:i+1])

  ;interpolate and trim to detector size
  ccd_flux  = interpol(smoothedflux,wave_new,wave_res,/SPLINE); linear interpolation
  ccd_wave  = wavelength
  ;read effective area and interpolate to cute wavelength
  length    = file_lines(file_eff)
  eff_area  = dblarr(5,length)
  openr,1,file_eff
  readf,1,eff_area
  close,1
  aeff        = interpol(eff_area[4,*],eff_area[0,*],wave_res,/SPLINE); linear interpolation
  ccd_count1  = dblarr(w_length/2)
  ;Effective area/QE
  ccd_count1  = ccd_flux*aeff*fwhm 
  ccd_count   = dblarr(w_length)
  noise       = dblarr(w_length)
  snr         = dblarr(w_length)   

  for i=0,w_length-1 do begin
    ccd_count[i]  = (ccd_count1[i/2]/2)*exptime ; assuming 2 resolution element 
    ;print,ccd_count1[i/2],ccd_count[i]
    noise[i]      = sqrt(ccd_count[i])+(width*(r_noise^2+(dark_noise*exptime)))
    snr[i]        = ccd_count[i]/noise[i]
  endfor


  ;user defined region
  userdefined=strsplit(infile.user_defined,' ',/EXTRACT)
  remove,[0,1],userdefined
  user_len=n_elements(userdefined)
  if evenodd(user_len) ne 0 then message,'Please recheck user defined inputs one value may be missing.'
  user_low=dblarr(user_len/2)
  user_high=dblarr(user_len/2)
  j=0
  for i=0,n_elements(userdefined)-2,2 do begin
    user_low[j]  = double(userdefined[i])
    user_high[j] = double(userdefined[i+1])
    j++
  endfor

  ;SNR calculations
  nx            = w_length
  spectra_1dwf  = ccd_count
  error_1dwf    = noise
  lenby3=fix(nx/3)
  lenby33=fix(nx-2*lenby3)
  st=0
  en=nx-1
  n_w=en-st+1
  x1=wavelength[st:en]
  y1=spectra_1dwf[st:en]
  dy1=error_1dwf[st:en]
  tp=trapz_error(x1,y1,dy1)
  snr_full=tp[0]/tp[1]
  st=0
  en=lenby3-1
  n_w=en-st+1
  x2=wavelength[st:en]
  y2=spectra_1dwf[st:en]
  dy2=error_1dwf[st:en]
  tp=trapz_error(x2,y2,dy2)
  snr_blue=tp[0]/tp[1]
  st=lenby3
  en=2*lenby3-1
  n_w=en-st+1
  x3=wavelength[st:en]
  y3=spectra_1dwf[st:en]
  dy3=error_1dwf[st:en]
  tp=trapz_error(x3,y3,dy3)
  snr_middle=tp[0]/tp[1]
  st=2*lenby3
  en=nx-1
  n_w=en-st+1
  x4=wavelength[st:en]
  y4=spectra_1dwf[st:en]
  dy4=error_1dwf[st:en]
  tp=trapz_error(x4,y4,dy4)
  snr_red=tp[0]/tp[1]
  ;near = Min(Abs(vector - number), index)
  ;2795,2802
  val=Min(Abs(wavelength - 2793), st)
  val2=Min(Abs(wavelength - 2805), en)
  n_w=en-st+1
  x5=wavelength[st:en]
  y5=spectra_1dwf[st:en]
  dy5=error_1dwf[st:en]
  tp=trapz_error(x5,y5,dy5)
  snr_mg2=tp[0]/tp[1]
  ;2852
  ;near = Min(Abs(vector - number), index)
  val=Min(Abs(wavelength - 2850), st)
  val2=Min(Abs(wavelength - 2854), en)
  n_w=en-st+1
  x6=wavelength[st:en]
  y6=spectra_1dwf[st:en]
  dy6=error_1dwf[st:en]
  tp=trapz_error(x6,y6,dy6)
  snr_mg1=tp[0]/tp[1]
  ;2585
  ;near = Min(Abs(vector - number), index)
  val=Min(Abs(wavelength - 2583), st)
  val2=Min(Abs(wavelength - 2587), en)
  n_w=en-st+1
  x7=wavelength[st:en]
  y7=spectra_1dwf[st:en]
  dy7=error_1dwf[st:en]
  tp=trapz_error(x7,y7,dy7)
  snr_fe2=tp[0]/tp[1]
  
  ;number of exposures in transit duration
  t_dur       = infile.transit_duration*3600. ;in seconds
  obs_time    = exptime+long(infile.readtime) ;in seconds
  obs_transit = t_dur/obs_time                ; number of observations in transit

  ;transit depth calculations per transit
  unc_tranit_full = 1./(sqrt(obs_transit)*snr_full)
  unc_tranit_blue = 1./(sqrt(obs_transit)*snr_blue)
  unc_tranit_middle = 1./(sqrt(obs_transit)*snr_middle)
  unc_tranit_red = 1./(sqrt(obs_transit)*snr_red)
  unc_tranit_mg1 = 1./(sqrt(obs_transit)*snr_mg1)
  unc_tranit_mg2 = 1./(sqrt(obs_transit)*snr_mg2)
  unc_tranit_fe2 = 1./(sqrt(obs_transit)*snr_fe2)

  ;transit depth calculations for n transit
  n=infile.number_of_transits
  n_unc_tranit_full   = unc_tranit_full/sqrt(n)
  n_unc_tranit_blue   = unc_tranit_blue/sqrt(n)
  n_unc_tranit_middle = unc_tranit_middle/sqrt(n)
  n_unc_tranit_red    = unc_tranit_red/sqrt(n)
  n_unc_tranit_mg1    = unc_tranit_mg1/sqrt(n)
  n_unc_tranit_mg2    = unc_tranit_mg2/sqrt(n)
  n_unc_tranit_fe2    = unc_tranit_fe2/sqrt(n)

  filename=infile.file_out+'cute_snr.txt'
  openw,1,filename
  printf,1,'-------------------------------------------------------------------------------'
  printf,1,'                      CUTE SNR CALCULATIONS'
  printf,1,'-------------------------------------------------------------------------------'
  printf,1,''
  printf,1,'Input Parameters'
  printf,1,''
  printf,1,'Stellar Parameters'
  printf,1,''
  printf,1,'Stellar Temperature(K)  : '+strtrim(string(infile.stellar_temperature),2)
  printf,1,'Stellar Radius (R_sun)  : '+strtrim(string(infile.stellar_radius),2)
  printf,1,'Stellar Magnitude (V)   : '+strtrim(string(infile.stellar_magnitude),2)
  printf,1,'Stellar Distance (mas)  : '+strtrim(string(infile.stellar_distance),2)
  printf,1,'Stellar Specrtral type  : '+strtrim(string(infile.spectral_type),2)
  printf,1,'Target RA (deg)         : '+strtrim(string(infile.Ra),2)
  printf,1,'Target Declination (deg): '+strtrim(string(infile.Dec),2)
  printf,1,''
  if infile.line_core_emission eq 1 then begin
    printf,1,'Line core emission added.'   
    printf,1,"logR'HK                 : "+strtrim(string(infile.logR),2)
  endif 
  if infile.add_ism_abs eq 1 then begin
    printf,1,'ISM absorption added for the following species'
    if (tag_exist(infile,'mg1_col')eq 1) then begin
      printf,1,'MgI column density      : '+strtrim(string(infile.mg1_col),2)
    endif else begin
      printf,1,'MgI column density (calculated).'
    endelse
    if (tag_exist(infile,'mg2_col')eq 1) then begin
      printf,1,'MgII column density     : '+strtrim(string(infile.mg2_col),2)
    endif else begin
      printf,1,'MgII column density (calculated).'
    endelse
    if (tag_exist(infile,'fe2_col')eq 1) then begin
      printf,1,'FeII column density     : '+strtrim(string(infile.fe2_col),2)
    endif else begin
      printf,1,'FeII column density (calculated).'
    endelse
  endif
  printf,1,''
  printf,1,'Instrument Parameters'
  printf,1,''
  printf,1,'Spectral Resolution (A) : '+strtrim(string(infile.spectral_resolution),2)
  printf,1,'Spectrum Width (pix)    : '+strtrim(string(infile.spectrum_width),2)
  printf,1,' in cross-disp direction  '
  printf,1,'Readout Noise (e/pix)   : '+strtrim(string(infile.readout_noise),2)
  printf,1,'Dark Noise (e/pix/s)    : '+strtrim(string(infile.dark_noise),2) 
  printf,1,'Exposure Time (s)       : '+strtrim(string(infile.exptime),2)
  printf,1,'Read Time (s)           : '+strtrim(string(infile.readtime),2)
  printf,1,''
  printf,1,'Input Files'
  printf,1,''
  in_file=infile.input_file 
  t=infile.stellar_temperature
  CASE StrUpCase(!Version.OS_Family) OF
    'WINDOWS': file=in_file+'\models\t'+String(t, Format='(I05)') +'g4.4\model.flx' ;assuming folder named by their temperature and file are named as model.flx ;WINDOWS
    'UNIX': file=in_file+'/models/t'+String(t, Format='(I05)') +'g4.4/model.flx' ;assuming folder named by their temperature and file are named as model.flx; UNIX.
  ENDCASE
  if file_test(file) ne 1 then t = t+100 ;above 8000K the steps is 200K
  CASE StrUpCase(!Version.OS_Family) OF
    'WINDOWS': file=in_file+'\models\t'+String(t, Format='(I05)') +'g4.4\model.flx' ;test again to see if temperature is not in range or is in steps of 100 or 200 ;WINDOWS
    'UNIX': file=in_file+'/models/t'+String(t, Format='(I05)') +'g4.4/model.flx' ;test again to see if temperature is not in range or is in steps of 100 or 200 UNIX.
  ENDCASE
  printf,1,'Stellar Model file      : '+strtrim(file,2)
  printf,1,'Wavelength file         : '+strtrim(string(infile.wave_file),2)
  printf,1,'Effective area file     : '+strtrim(string(infile.file_eff),2)
  printf,1,'-------------------------------------------------------------------------------'
  printf,1,'                        CUTE SNR OUTPUTS'
  printf,1,'-------------------------------------------------------------------------------'
  printf,1,'Wavelength Region [A]           SNR          Uncertanity in Transit Depth [ppm]'
  printf,1,'                                             1 Transit    '+strtrim(n,2)+' Transits'  
  printf,1,'-------------------------------------------------------------------------------'  
  printf,1,'Full Band['+strtrim(string(x1[0],Format='(F7.2)'),2)+'-'+strtrim(string(x1[n_elements(x1)-1],Format='(F7.2)'),2)+']      '$
        +strtrim(string(snr_full),2)+'    '+strtrim(string(unc_tranit_full*1E6),2)+'    '+strtrim(string(n_unc_tranit_full*1E6),2)
  printf,1,'Lower Band['+strtrim(string(x2[0],Format='(F7.2)'),2)+'-'+strtrim(string(x2[n_elements(x2)-1],Format='(F7.2)'),2)+']     '$
        +strtrim(string(snr_blue),2)+'    '+strtrim(string(unc_tranit_blue*1E6),2)+'    '+strtrim(string(n_unc_tranit_blue*1E6),2)
  printf,1,'Middle Band['+strtrim(string(x3[0],Format='(F7.2)'),2)+'-'+strtrim(string(x3[n_elements(x3)-1],Format='(F7.2)'),2)+']    '$
        +strtrim(string(snr_middle),2)+'    '+strtrim(string(unc_tranit_middle*1E6),2)+'    '+strtrim(string(n_unc_tranit_middle*1E6),2)
  printf,1,'Upper Band['+strtrim(string(x4[0],Format='(F7.2)'),2)+'-'+strtrim(string(x4[n_elements(x4)-1],Format='(F7.2)'),2)+']     '$
        +strtrim(string(snr_red),2)+'    '+strtrim(string(unc_tranit_red*1E6),2)+'    '+strtrim(string(n_unc_tranit_red*1E6),2)
  printf,1,'MgII Band['+strtrim(string(x5[0],Format='(F7.2)'),2)+'-'+strtrim(string(x5[n_elements(x5)-1],Format='(F7.2)'),2)+']      '$
        +strtrim(string(snr_mg2),2)+'    '+strtrim(string(unc_tranit_mg2*1E6),2)+'    '+strtrim(string(n_unc_tranit_mg2*1E6),2)
  printf,1,'MgI Band['+strtrim(string(x6[0],Format='(F7.2)'),2)+'-'+strtrim(string(x6[n_elements(x6)-1],Format='(F7.2)'),2)+']       '$
        +strtrim(string(snr_mg1),2)+'    '+strtrim(string(unc_tranit_mg1*1E6),2)+'    '+strtrim(string(n_unc_tranit_mg1*1E6),2)
  printf,1,'FeII Band['+strtrim(string(x7[0],Format='(F7.2)'),2)+'-'+strtrim(string(x7[n_elements(x7)-1],Format='(F7.2)'),2)+']      '$
        +strtrim(string(snr_fe2),2)+'    '+strtrim(string(unc_tranit_fe2*1E6),2)+'    '+strtrim(string(n_unc_tranit_fe2*1E6),2) 
  for i=0,(user_len/2)-1 do begin
    ;near = Min(Abs(vector - number), index)
    ud_low    = user_low[i]
    ud_high   = user_high[i]
    val       = Min(Abs(wavelength - ud_low), st)
    val2      = Min(Abs(wavelength - ud_high), en)
    n_w       = en-st+1
    x_ud   = wavelength[st:en]
    y_ud   = spectra_1dwf[st:en]
    dy_ud  = error_1dwf[st:en]
    tp        = trapz_error(x_ud,y_ud,dy_ud)
    snr_ud = tp[0]/tp[1]
    unc_tranit_ud    = 1./(sqrt(obs_transit)*snr_ud)
    n_unc_tranit_ud   = unc_tranit_ud/sqrt(n)
    printf,1,'User Band '+strtrim(string(i),2)+'['+strtrim(string(x_ud[0],Format='(F7.2)'),2)+'-'+strtrim(string(x_ud[n_elements(x_ud)-1],Format='(F7.2)'),2)+']    '$
         +strtrim(string(snr_ud),2)+'    '+strtrim(string(unc_tranit_ud*1E6),2)+'    '+strtrim(string(n_unc_tranit_ud*1E6),2)  
  endfor                
  printf,1,'-------------------------------------------------------------------------------'
  printf,1,'Full Spectrum SNR'
  printf,1,''
  printf,1,'Wavelength [A]  Flux [counts]      Noise[counts]      SNR'
  for i=0,w_length-1 do printf,1,FORMAT='(F12.2,1x,3(F17.9,1x))',wavelength[i],ccd_count[i],noise[i],snr[i]
  close,1            

  ;writecol,infile.file_out+'cute_snr.txt',wavelength[i],ccd_count[i],noise[i],snr[i],fmt='(4(F17.9,1x))
  ;pson,filename=infile.file_out+'cute_snr_plot.ps'
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
end

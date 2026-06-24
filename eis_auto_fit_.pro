

PRO eis_auto_fit_, windata, fitdata, template=template, wvlpix=wvlpix, $
                  xrange=xrange, yrange=yrange, quiet=quiet, $
                  offset=offset, uniform_backg=uniform_backg, $
                  wvl_select=wvl_select,refwvl=refwvl,perpixel=perpixel, $
                  iexp=iexp, mask=mask, FUNCTION_name=FUNCTION_name

;+
; NAME
;
;   EIS_AUTO_FIT
;
; PROJECT
;
;   Hinode/EIS
;
; EXPLANATION
;
;   This routine automatically fits single or multiple Gaussians to
;   two-dimensional spatial arrays of Hinode/EIS spectra. It
;   supercedes the earlier routines eis_auto_fit and
;   eis_auto_fit_gen. 
;
;   eis_auto_fit is designed to handle modified windata structures
;   produced by, e.g., eis_bin_windata and eis_join_windata, in
;   addition to the standard structures.
;
;   Please consult the online tutorial (available on the EIS wiki) for
;   help with using eis_auto_fit.
;
;   For both multi-Gauss and single Gauss fits, the code makes use of
;   the MPFIT parinfo structure to constrain the range of allowed fit
;   parameters. For example a line's amplitude is forced to be
;   > 0. If one of these limits is reached, then EIS_AUTO_FIT
;   forces the fit parameters for that pixel to be set to zero and a
;   warning message is printed. If a lot of cases are found then the
;   user should think about manually changing these limits in the
;   TEMPLATE structure.
;
; INPUTS
;
;   WINDATA   A window data structure in the format created by
;             EIS_GETWINDATA. 
;
; OPTIONAL INPUTS
;
;   TEMPLATE  A structure that specifies a template for the
;             Gaussian fit - see the routine EIS_MAKE_FIT_TEMPLATE
;             for more details. Note that it is essential to specify
;             this structure if you're doing a multi-Gauss fit,
;             but it is optional for single Gauss fits.
;
;   WVL_SELECT A structure created by EIS_WVL_SELECT that specifies
;              which wavelengths are to be included in the
;              fitting. See EIS_WVL_SELECT for more details.
;
;   WVLPIX   By default eis_auto_fit_new assumes the entire wavelength range in
;            the data array will be used for fitting. The range can be reduced
;            by specifying the 2 element array WVLPIX. E.g., setting
;            WVLPIX=[194.9,195.4] will use only the wavlength range 194.9
;            to 195.4 angstroms (if the Fe XII 195 line is being
;            fit). If the OFFSET array is being input to
;            EIS_AUTO_FIT_NEW, then it is *highly recommended* that
;            you do not use WVLPIX, but use WVL_SELECT instead.
;
;   OFFSET   A 2D array of same size as the spatial dimensions of
;            WINDATA that gives the wavelength offset that applies at
;            each pixel due to the EIS slit tilt and orbit variation.
;
;   XRANGE   Allows a sub-range of the X-values to be fitted. E.g., setting
;            equal to [10,30] means that only X-values from 10 to 30
;            (inclusive) are fit.
;
;   YRANGE   Allows a sub-range of the Y-values to be fitted.  E.g., setting
;            equal to [10,30] means that only Y-values from 10 to 30
;            (inclusive) are fit.
;
;     IEXP   This is only active if the data-set has nexp_prp>1
;            (i.e., more than one exposure per raster
;            position). It is an index that selects which exposure
;            to use. By default the routine selects the exposure
;            with the longest exposure time.
;
;     MASK   A 2D byte array of the same size as the spatial
;            dimensions of WINDATA. If the value of a pixel is 0, then
;            that spatial pixel will be excluded from the fit, whereas
;            if it's 1 it will be included.
;
;     Function_Name: Specify the name of a function for fitting the lines. The
;               function must be defined by three parameters that correspond
;               to peak, centroid and width. If not specified, then the
;               function gauss_sg is used.
;
; KEYWORD PARAMETERS:
;
;   UNIFORM_BACKG By default, the routine performs a linear fit to the
;                 spectrum background. By setting this keyword the background
;                 is assumed flat.
;
;   QUIET    If set, then do not print any text information to the IDL
;            window.
;
;   PERPIXEL   EIS_AUTO_FIT assumes that the intensity is given in
;              "per-Angstrom" units. By setting /PERPIXEL "per-pixel"
;              units are assumed instead. 
;
; OUTPUTS:
;
;   FITDATA  A structure containing the fits to the emission
;            line(s). See the section 'PROGRAMMING NOTES' for the list
;            of tags. To extract information such as velocities and
;            line widths from FITDATA, it is recommended that the
;            routine EIS_GET_FITDATA is used.
;
; PROGRAMMING NOTES
;
;   The output structure 'FITDATA' has the following tags:
;
;   .aa  The fit parameters for each spatial pixel. It has the
;        dimensions (nparams,nx,ny). The line peaks (amplitudes) are
;        stored in positions 0, 3, 6, etc., the centroids in positions
;        1, 4, 7, etc., and the line widths (Gaussian widths, not
;        FWHM) in 2, 5, 8, etc. The parameters for the background
;        (either 1, if /uniform_backg set, or 2) are appended to the line
;        parameters. E.g., if there are two lines, then the background
;        parameters will be in positions 6 and 7. If a fit is not
;        performed for a spatial pixel, then the AA values for that
;        pixel will be set to 0 (not MISSING). See the tag 'bad_pix'
;        for more details about pixels with bad fits.
;
;   .sigmaa  The 1-sigma errors on the parameters stored in .aa. If a
;            fit has not been performed for a pixel then sigmaa is set
;            to 0. See the tag 'bad_pix' for more details about pixels
;            with bad fits.
;
;   .int The intensity arrays for each line. Has dimensions
;        (nlines,nx,ny). The intensity is computed from the line peaks
;        and widths stored in the .aa tag.
;
;   .interr The 1-sigma error on the intensities stored in .int. Has
;           dimensions (nlines,nx,ny).
;
;   .chi2  Contains the reduced chi^2 values from the fits. Has
;          dimensions (nx,ny)
;
;   .nx    Size of the fit arrays in the Solar-X dimension.
;
;   .wd_nx The original size of the raster arrays of the observation
;          (this can be different to .nx)
;
;   .ny    Size of the fit arrays in the Solar-Y dimension.
;
;   .wd_ny The original size of the raster arrays of the observation
;          (this can be different to .ny).
;
;   .offset The array of wavelength offsets generated by, e.g.,
;           EIS_WAVE_CORR. If not specified, then set to zero. Has
;           dimensions (nx,ny) 
;
;   .refwvl Contains reference wavelengths for each line of fit. Has
;           dimension nlines.
;
;   .x_bg1  The lower wavelength at which the background is
;           specified. See below for more details about the
;           background. Has dimensions (nx,ny).
;
;   .x_bg2  The upper wavelength at which the background is
;           specified. See below for more details about the
;           background. Has dimensions (nx,ny).
;
;   .exp_start_times  The start time of each exposure. Has dimensions
;                     nx. 
;
;   .missing The value of missing data in the array (inherited from
;            WINDATA). 
;
;   .xrange The X-range (in pixel indices) used for the fit. E.g.,
;           [40,59] implies pixels 40 to 59 from within the specified
;           WINDATA arrays were used. (See XRANGE optional input.)
;
;   .yrange The Y-range (in pixel indices) used for the fit. E.g.,
;           [40,59] implies pixels 40 to 59 from within the specified
;           WINDATA arrays were used. (See YRANGE optional input.)
;
;   .ngauss The number of Gaussians (i.e., emission lines) used in the
;           fit. 
;
;   .back   The number of background parameters used in the fit.
;
;   .slit_ind The EIS slit index for the observation (0 - 1", 2 -
;             2").
;
;   .ybin   The size of Y-bins in the WINDATA structure. Usually this
;           is 1, but see the routine eis_bin_windata for more
;           details. 
;
;   .yip    The Y-pixel on the detector corresponding to the bottom
;           pixel of the raster.
;
;   .raster Takes the value 1 or 0 for if the observation is a raster
;           or sit-and-stare, respectively.
;
;   .start_pix A two element array specifying the start pixel used by
;              the routine EIS_BIN_WINDATA.
;
;   .date_obs  Start time of observation (obtained from WINDATA
;              header). 
;
;   .bad_pix  This is a byte array of same size as INT that indicates
;             if there was something wrong with the spatial pixel that
;             prevented a fit from being performed. A value
;             of 1 is because the spectrum is too noisy. A value of 2
;             is because the number of pixels 
;             being fit is less than the number of free parameters. A
;             special case of the latter is when the whole spectrum is
;             missing (bad_pix=3). A value of 4 indicates the offset
;             value for that pixel was missing. A value of 5 indicates
;             that the fit failed (sigmaa ended up undefined).
;
;    .time_stamp  String that gives time at which file was created.
;
;    .iexp_prp If nexp_prp>1 then iexp_prp gives exposure index (between 0
;              and nexp_prp-1. If nexp_prp=1, then iexp_prp=-1.
;
;    .instrume    Name of instrument to which data apply ('EIS' or
;                 'CDS', at present).
;
;    .function String giving name of fit function. E.g., 'gauss_sg'.
;
;  
;   The background in the fit is defined to be a straight line with
;   values specified at the lowest and highest wavelengths of the fit
;   range. I.e., suppose the region from 199 to 200 A is being fitted,
;   then the straight line is defined by the values at 199 and 200
;   A. This becomes a little more complicated by the wavelength
;   offsets stored in the OFFSET input, so the wavelength values used
;   for each spatial pixel are stored in the output tags .x_bg1 and
;   .x_bg2. The actual fit values at these wavelength are stored in
;   the fit parameter tag, .aa.
;
; CALLS
;
;    EIS_TEMPLATE_PARINFO, EIS_READ_TEMPLATE, EIS_FIT_FUNCTION
;
; HISTORY
;
;    Version 1, 2-Feb-2010, Peter Young
;      code adapted from eis_auto_fit and eis_auto_fit_gen
;    Version 2, 15-Apr-2010, Peter Young
;      modified use of bad_pix (see programming notes); routine no
;      longer complains about parameter limits when a fit has not been
;      performed; missing data in the .int, .err and .chi2 tags are
;      now correctly set to windata.missing; check on noisy data has
;      been modified: max(spec)-median(spec) must now be >
;      2*max(error); also now prints statistics on pixels for which
;      fits were not possible.
;    Version 3, 23-Apr-2010, Peter Young
;      If the offset value at a pixel is missing, then a fit is not
;      attempted (and bad_pix set to 4); modified check on offset in
;      case it is a single element.
;    Version 4, 12-Apr-2011, Peter Young
;      The REFWVL optional input has been disabled. If OFFSET is not
;      specifed, then the routine now uses WINDATA.WAVE_CORR.
;    Version 5, 3-May-2011, Peter Young
;      Initial parameters for template case now correctly scaled
;      (search for "3-May-11" in code).
;    Version 6, 9-May-2011, Peter Young
;      Added some checks on template to make sure that there is at
;      least one free parameter for each of line peak, line width and
;      centroid. 
;    Version 7, 11-May-2011, Peter Young
;      Added scale, xcen and ycen tags to fitdata output.
;    Version 8, 8-Apr-2012, Peter Young
;      Now warns about parameter limits being hit.
;    Version 9, 26-Apr-2012, Peter Young
;      Added new bad_pix=5 case when fit fails to be performed; added
;      check on template.lines.wid_lim as this was a common source of
;      error. 
;    Version 10, 12-Jul-2012, Peter Young
;      Added check on template to make sure it's consistent
;      with windata. 
;    Version 11, 15-Aug-2012, Peter Young
;      A few minor changes to text and print statements; main code is
;      unchanged. 
;    Version 12, 8-Jan-2013, Peter Young
;      Added time_stamp to FITDATA.
;    Version 13, 4-Apr-2013, Peter Young
;      I'm attempting to make the auto_fit suite of routines
;      work with CDS data, so I've added the INSTRUME tag to
;      the output.
;    Version 14, 23-Jun-2014, Peter Young
;      Added /perpixel option in order for routine to work with IRIS
;      data. 
;    Version 15, 18-Aug-2014, Peter Young
;      Modified the 'instrume' tag so that it is either CDS, EIS or
;      IRIS; fixed text format problem for IRIS; switched to using the
;      'progress' routine; modified initial width values for different
;      instruments. 
;     Version 16, 2-Mar-2015, Peter Young
;      Fixed bug for nexp_prp>1 data; added iexp_prp tag to output structure.
;     Version 17, 24-Jun-2015, Peter Young
;      Modified check on strength of signal compared to error (search
;      for 24-Jun-2015 in code).
;     Version 18, 13-Nov-2015, Peter Young
;      The /quiet keyword wasn't stopping all of the print output, so
;      I've fixed this now.
;     Version 19, 5-Feb-2018, Peter Young
;      Added the MASK optional input.
;     Version 20, 13-Feb-2018, Peter Young
;      Fixed bug when /perpixel was used. I was dividing by
;      wvl_factor, which made the error too large.
;     Version 21, 10-May-2020, Peter Young
;      Updated my contact details; no change to code.
;     Version 22, 24-Aug-2020, Peter Young
;      Fixed bug for when wvl_select is not specified; only a problem
;      if wave_corr was a 9x9 array or smaller (so very rarely).
;     Version 23, 18-Dec-2020, Peter Young
;      Caught bug when all of the spectral pixels are missing for one
;      of the spatial pixels.
;     Version 24, 24-Jan-2022, Peter Young
;      Added a fix for SPICE, but routine still does not work with
;      SPICE data (errors are not currently defined for SPICE, and
;      need to handle NaN values for missing).
;     Version 25, 05-Nov-2024, Peter Young
;      Added optional input function_name=; added function tag to output
;      structure.
;     Version 26, 14-May-2025, Peter Young
;      If the intensity units for EIS are DN or photons, then
;      automatically set the perpixel keyword. I removed the warning
;      about the perpixel keyword (see v.20) since it's quite old now.
;     Version 27, 18-Feb-2026, Peter Young
;      Introduced nfree for the number of free parameters, computed
;      using the template data. This is used to determine if there are
;      enough data points for the fit to work.
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use: IDL> eis_auto_fit, windata, fitdata [, template=, wvlpix=, '
  print,'                     xrange=, yrange=, offset=, wvl_select= '
  print,'                     /quiet, /uniform_backg, mask= ]'
  return
ENDIF

IF NOT keyword_set(quiet) THEN BEGIN
  print,''
  print,' EIS_AUTO_FIT was written by Peter Young (GSFC).'
  print,' Please report any errors to peter.r.young@nasa.gov.'
  print,''
ENDIF 


t1=systime(1)

IF n_tags(template) NE 0 AND keyword_set(uniform_backg) THEN BEGIN
  print,'%EIS_AUTO_FIT: instead of using the /UNIFORM_BACKG, when giving a fit template the user should set:'
  print,'                 IDL> template.nback=1'
  print,'               Please perform this step and then call eis_auto_fit without the /UNIFORM_BACKG keyword.'
  print,''
  return
ENDIF 

;
; The intensity units of EIS data are either erg/cm2/s/sr/A, DN, or Photons.
; In the latter two cases the perpixel keyword needs to be set.
;
IF tag_exist(windata,'units') AND windata.hdr.instrume EQ 'EIS' THEN BEGIN
  IF n_elements(perpixel) EQ 0 THEN BEGIN 
    IF windata.units EQ 'DN' THEN perpixel=1b
    IF windata.units EQ 'Photon-Events' THEN perpixel=1b
  ENDIF 
ENDIF 

nl=windata.nl
nx=windata.nx
ny=windata.ny
wvl=windata.wvl
nexp_prp=windata.hdr.nexp_prp
IF nexp_prp GT 1 THEN BEGIN
 ;
  exp_times=reform(windata.exposure_time[0,*])
  nexp=n_elements(exp_times)
  IF n_elements(iexp) EQ 0 THEN BEGIN
    getmax=max(exp_times,iexp)
    print,'% EIS_AUTO_FIT: this data-set has nexp_prp>1. The exposure times are:'
    FOR i=0,nexp-1 DO BEGIN
      print,format='(i10,f8.1," s")',i,exp_times[i]
    ENDFOR 
    print,' The largest exposure time will be used. Use iexp= to select a different exposure index.'
    print,''
  ENDIF ELSE BEGIN
   ;
    IF iexp LT 0 OR iexp GE nexp_prp THEN BEGIN
      print,'% EIS_AUTO_FIT: iexp should take a value between 0 and '+trim(nexp_prp)+'.'
      return
    ENDIF 
   ;
  ENDELSE 
  int=windata.int[*,*,*,iexp]
  err=windata.err[*,*,*,iexp]
  exp_start_times=windata.time[*,iexp]
ENDIF ELSE BEGIN 
  int=windata.int
  err=windata.err
  exp_start_times=windata.time
  iexp=-1
ENDELSE 
miss=windata.missing


;
; The following implements the MASK input.
;  Any spatial pixels flagged with 0 in MASK will set the
;  corresponding spatial pixels in INT and ERR to be missing.
;
IF n_elements(mask) NE 0 THEN BEGIN
  s=size(mask,/dim)
  test1=s[0] NE nx OR s[1] NE ny
  test2= min(mask) LT 0 OR max(mask) GT 1 
  IF test1 THEN BEGIN
    print,'% EIS_AUTO_FIT: the dimensions of MASK do not match the spatial dimensions of WINDATA.'
    print,format='("                X:",i4,",",i4," Y:",i4,",",i4)',s[0],nx,s[1],ny
    print,'                MASK will be ignored.'
  ENDIF
  IF test2 THEN BEGIN
    print,"% EIS_AUTO_FIT: MASK must contain only 0's or 1's. The input will be ignored."
  ENDIF
  IF test1 EQ 0 AND test2 EQ 0 THEN BEGIN
    mask_arr=bytarr(nl,nx,ny)
    FOR i=0,nl-1 DO mask_arr[i,*,*]=mask
    k=where(mask_arr EQ 0,nk)
    IF nk NE 0 THEN BEGIN
      int[k]=miss
      err[k]=miss
    ENDIF 
  ENDIF 
ENDIF


;
; Check to make sure TEMPLATE is consistent with WINDATA
; 
IF n_tags(template) NE 0 THEN BEGIN 
  k=where(template.lines.centroid GE min(windata.wvl) AND template.lines.centroid LE max(windata.wvl),nk)
  IF nk EQ 0 THEN BEGIN
    print,'%EIS_AUTO_FIT:  The wavelengths in TEMPLATE are not consistent with WINDATA. Please check '
    print,'                your inputs. Returning...'
    return
  ENDIF
ENDIF 


;
; The input TEMPLATE can either be the structure returned by
; eis_read_template, or it can be the name of a fit template file
; (which will be read by eis_read_fit_template).
;
IF n_elements(template) NE 0 THEN BEGIN
  IF n_tags(template) EQ 0 THEN BEGIN
    chck=file_exist(template)
    IF chck EQ 1 THEN BEGIN
      template=eis_read_template(template)
    ENDIF
  ENDIF
 ;
 ; Perform some sanity checks on TEMPLATE
 ;
  k=where(template.lines.cen_tie EQ -1,nk)
  nlines=template.ngauss
  IF nk EQ 0 THEN BEGIN
    print,'%EIS_AUTO_FIT: Problem with the TEMPLATE input!'
    print,'   There are no free centroid parameters in your fit template:'
    print,'      Line  cen_tie value'
    FOR i=0,nlines-1 DO print,format='(2i10)',i,template.lines[i].cen_tie
    print,'   At least one cen_tie value should be -1. Please check or re-generate your TEMPLATE.'
    print,'   EIS Software Note #16 has more details on how TEMPLATE is specified.'
    print,'   Returning...'
    fit=-1
    return
  ENDIF 
 ;
  k=where(template.lines.wid_tie EQ -1,nk)
  nlines=template.ngauss
  IF nk EQ 0 THEN BEGIN
    print,'%EIS_AUTO_FIT: Problem with the TEMPLATE input!'
    print,'   There are no free line width parameters in your fit template:'
    print,'      Line  wid_tie value'
    FOR i=0,nlines-1 DO print,format='(2i10)',i,template.lines[i].wid_tie
    print,'   At least one wid_tie value should be -1. Please check or re-generate your TEMPLATE.'
    print,'   EIS Software Note #16 has more details on how TEMPLATE is specified.'
    print,'   Returning...'
    fit=-1
    return
  ENDIF 
 ;
  k=where(template.lines.peak_tie EQ -1,nk)
  nlines=template.ngauss
  IF nk EQ 0 THEN BEGIN
    print,'%EIS_AUTO_FIT: Problem with the TEMPLATE input!'
    print,'   There are no free line peak parameters in your fit template:'
    print,'      Line  peak_tie value'
    FOR i=0,nlines-1 DO print,format='(2i10)',i,template.lines[i].peak_tie
    print,'   At least one wid_tie value should be -1. Please check or re-generate your TEMPLATE.'
    print,'   EIS Software Note #16 has more details on how TEMPLATE is specified.'
    print,'   Returning...'
    fit=-1
    return
  ENDIF 
 ;
 ; This is a common problem (at least for me!)...accidentally setting wid_lim to a
 ; single value rather than fltarr(2).
 ;
  FOR i=0,nlines-1 DO BEGIN
    wid_lim=template.lines[i].wid_lim
    IF wid_lim[0] EQ wid_lim[1] AND wid_lim[0] NE template.null_value THEN BEGIN
      print,'%EIS_AUTO_FIT: Problem with template.lines['+trim(i)+'].wid_lim !'
      print,'               Values are identical: ',wid_lim
      print,'               Please specify a width range.  Returning...'
      fit=-1
      return
    ENDIF
  ENDFOR 
ENDIF


;
; If TEMPLATE exists, then use this to determine number of Gaussians. 
; If it doesn't exist, then assume a single Gaussian is being fit.
;
; nparams is the total number of parameters for Gaussians plus background.
; nfree is the number of free parameters, and is computed using the tied
;       parameter information in template.
;
IF n_tags(template) NE 0 THEN BEGIN
  ngauss=n_elements(template.lines)
  IF keyword_set(uniform_backg) THEN BEGIN
    nback=1 
  ENDIF ELSE BEGIN
    nback=template.nback
  ENDELSE
  ;
  ; Only if the "tie" value is -1 is the parameter free.
  k=where(template.lines.peak_tie EQ -1,nparams_peak)
  k=where(template.lines.cen_tie EQ -1,nparams_cen)
  k=where(template.lines.wid_tie EQ -1,nparams_wid)
  ;
  nfree=nparams_peak+nparams_cen+nparams_wid+nback
  nparams=3*ngauss+nback
 ;
 ; Obtain parinfo structure from TEMPLATE
 ;
  parinfo=eis_template_parinfo(template)
ENDIF ELSE BEGIN
  ngauss=1
  IF keyword_set(uniform_backg) THEN BEGIN
    nback=1
    nparams=4
  ENDIF ELSE BEGIN
    nback=2
    nparams=5
  ENDELSE
ENDELSE 

;
; Using parameters defined above, create the fit function.
;
IF n_elements(template) NE 0 THEN BEGIN
  expr=eis_fit_FUNCTION(template, FUNCTION_name=FUNCTION_name)
 ;
 ; The following checks for those parameters that are tied, and sets
 ; them to be fixed quantities.
 ; Note that n_free gives the number of free parameters
 ;
  n_free=nparams
  FOR i=0,ngauss-1 DO BEGIN
    peak_tie=template.lines[i].peak_tie
    cen_tie=template.lines[i].cen_tie
    wid_tie=template.lines[i].wid_tie
    IF (peak_tie GE 0) AND (peak_tie NE i) THEN BEGIN
      parinfo[i*3].fixed=1
      n_free=n_free-1
    ENDIF 
    IF (cen_tie GE 0) AND (cen_tie NE i) THEN BEGIN
      parinfo[i*3+1].fixed=1
      n_free=n_free-1
    ENDIF 
    IF (wid_tie GE 0) AND (wid_tie NE i) THEN BEGIN
      parinfo[i*3+2].fixed=1
      n_free=n_free-1
    ENDIF 
  ENDFOR
ENDIF ELSE BEGIN
  expr=''
  IF n_elements(FUNCTION_name) EQ 0 THEN func='gauss_sg' ELSE func=FUNCTION_name
  FOR i=0,ngauss-1 DO BEGIN
    expr=expr+func+'(x,p['+trim(i*3)+':'+trim(i*3+2)+']) + '
  ENDFOR 
  IF nback EQ 2 THEN BEGIN
    expr=expr+'line_sg(x,p['+trim(ngauss*3)+':'+trim(ngauss*3+1)+'])'
  ENDIF ELSE BEGIN
    expr=expr+'p['+trim(ngauss*3)+']'
  ENDELSE
  n_free=ngauss*3+nback
ENDELSE 

;
; Set any missing pixels (flagged in the error array) to be missing in
; the intensity array.
;
i=where(err EQ miss,ni)
IF ni NE 0 THEN int[i]=miss

nwvl=n_elements(wvl)


IF n_elements(xrange) NE 0 THEN BEGIN
  nx=xrange[1]-xrange[0]+1
  exp_start_times=exp_start_times[xrange[0]:xrange[1]]
 ;
  scale=windata.scale
  origx=windata.xcen - windata.nx/2.*scale[0]
  xcen=origx + xrange[0]*scale[0] + (xrange[1]-xrange[0]+1)/2.*scale[0]
ENDIF ELSE BEGIN
  xrange=[0,nx-1]
  xcen=windata.xcen
ENDELSE

IF n_elements(yrange) NE 0 THEN BEGIN
  ny=yrange[1]-yrange[0]+1 
 ;
  scale=windata.scale
  origy=windata.ycen - windata.ny/2.*scale[1]
  ycen=origy + yrange[0]*scale[1] + (yrange[1]-yrange[0]+1)/2.*scale[1]
ENDIF ELSE BEGIN
  yrange=[0,ny-1]
  ycen=windata.ycen
ENDELSE 


;
; This sets up the array that will contain the reference wavelengths
; (which will be defined after the fits have been done).
;
refwvl=dblarr(ngauss)



IF NOT tag_exist(windata,'wave_corr') AND n_elements(offset) EQ 0 THEN BEGIN
  print,'%EIS_AUTO_FIT: Your WINDATA structure is out-of-date and does not contain the WAVE_CORR tag.'
  print,'Please re-generate the structure with EIS_GETWINDATA before calling EIS_AUTO_FIT.'
  print,'Returning...'
  return
ENDIF 

IF n_elements(offset) EQ 0 THEN BEGIN
  offset=windata.wave_corr
ENDIF ELSE BEGIN
  siz=size(offset)
  IF siz[1] NE windata.nx OR siz[2] NE windata.ny THEN BEGIN
    print,'%EIS_AUTO_FIT: the input OFFSET= has been specified, but the dimensions do not match those of WINDATA.'
    print,'Please check your data.  Returning...'
    return
  ENDIF 
ENDELSE 

;
; offset_use is the offset array used in the code; it can be a
; sub-array of the input offset if xrange and/or yrange are set.
;
offset_use=offset[xrange[0]:xrange[1],yrange[0]:yrange[1]]

;
; 'raster' indicates whether the data-set represents a raster
; (raster=1)  or sit-and-stare (raster=0). This is needed when the
; orbit correction is done by eis_update_fitdata.pro.
;
IF windata.hdr.nraster GT 1 THEN raster=1 ELSE raster=0

;
; 'start_pix' is a tag introduced by eis_bin_windata and is
; non-standard therefore first need to check that tag exists in
; structure
;
IF tag_exist(windata,'start_pix') EQ 1 THEN BEGIN
  start_pix=windata.start_pix
ENDIF ELSE BEGIN
  start_pix=[0,0]
ENDELSE

;
; If WVL_SELECT was not specified then create the structure, using
; WVL_PIX if this has been specified.
;
IF n_tags(wvl_select) EQ 0 THEN BEGIN
  wvl_select={n: 1, min: fltarr(10), max: fltarr(10)}
  IF n_elements(wvlpix) NE 0 THEN BEGIN
    wvl_select.min[0]=wvlpix[0]
    wvl_select.max[0]=wvlpix[1]
  ENDIF ELSE BEGIN
    wvl_select.min[0]=min(windata.wvl)-max(offset_use)
    wvl_select.max[0]=max(windata.wvl)-min(offset_use)
  ENDELSE
ENDIF 


;
; This extracts the name of the instrument. It should be 'EIS', 'CDS'
; or 'IRIS'. 
;
IF tag_exist(windata.hdr,'telescop') THEN BEGIN
  telescop=windata.hdr.telescop
  IF trim(strlowcase(telescop)) EQ 'iris' THEN BEGIN
    instrume='IRIS'
  ENDIF ELSE BEGIN 
    instrume=trim(windata.hdr.instrume)
  ENDELSE 
ENDIF


;
; Define the FITDATA structure
;
IF n_elements(FUNCTION_name) EQ 0 THEN func='gauss_sg' ELSE func=FUNCTION_name
fitdata={exp_start_times: exp_start_times, $
         missing: miss, $
         xrange: xrange, $
         yrange: yrange, $
         ngauss: ngauss, $
         nback: nback, $
         refwvl: refwvl, $
         nx: nx, $
         ny: ny, $
         offset: offset_use, $
         aa: fltarr(nparams,nx,ny), $
         sigmaa: fltarr(nparams,nx,ny), $
         int: fltarr(ngauss,nx,ny)+miss, $
         interr: fltarr(ngauss,nx,ny)+miss, $
         chi2: fltarr(nx,ny)+miss, $
         x_bg1: fltarr(nx,ny), $
         x_bg2: fltarr(nx,ny), $
         slit_ind: windata.hdr.slit_ind, $
         ybin: windata.scale[1], $
         wd_nx: windata.hdr.nexp, $
         wd_ny: windata.hdr.yw, $
         raster: raster, $
         start_pix: start_pix, $
         yip: windata.hdr.yws, $
         bad_pix: bytarr(nx,ny), $
         wvl_select: wvl_select, $
         date_obs: windata.hdr.date_obs, $
         xcen: xcen, $
         ycen: ycen, $
         scale: windata.scale, $
         time_stamp: '', $
         iexp_prp: iexp, $
         FUNCTION_name: func, $
         instrume: instrume }


;
; If template hasn't been specified, then I need to set initial
; parameters for the line width, which varies significantly between instruments.
;
CASE instrume OF
  'IRIS': BEGIN
    init_wid=0.03
    wid_range=[0.005,0.55]
  END 
  'EIS': BEGIN
    init_wid=0.03
    wid_range=[0.02,0.10]
  END 
  'CDS': BEGIN 
    init_wid=0.2
    wid_range=[0.1,0.5]
  END
  'SPICE': BEGIN
    init_wid=0.9
    wid_range=[0.4,2.0]
  END 
ENDCASE 

IF NOT keyword_set(quiet) THEN progress,0.,/reset

;
; The two for loops below go through each spatial pixel and perform
; the fitting.
;
FOR i=xrange[0],xrange[1] DO BEGIN
  FOR k=yrange[0],yrange[1] DO BEGIN

    i0=i-xrange[0]
    k0=k-yrange[0]

    offset_val=offset_use[i0,k0]

    IF offset_val NE windata.missing THEN BEGIN 

      xx=wvl - offset_val

    IF n_tags(wvl_select) NE 0 THEN BEGIN
      pix_all=bytarr(nwvl)
      nws=wvl_select.n
      FOR iw=0,nws-1 DO BEGIN
        chck=where(xx GE wvl_select.min[iw] AND xx LE wvl_select.max[iw],nchck)
        IF nchck NE 0 THEN pix_all[chck]=1b
      ENDFOR 
    ENDIF 

   ;
   ; PRY, 18-Dec-2020
   ; SWTCH is used here to determine whether to go ahead with the
   ; fitting.
   ;
    swtch=1b
    ipix=where(pix_all EQ 1,n_ipix)
    IF n_ipix EQ 0 THEN BEGIN
       swtch=0b
    ENDIF ELSE BEGIN 
       xx=xx[ipix]
       yy=int[ipix,i,k]
       ee=err[ipix,i,k]
      ;
       IF (max(yy)-median(yy)) LE 2.0*median(ee) THEN swtch=0b
    ENDELSE 
   
   ;
   ; The check below is to stop the routine wasting time trying to fit
   ; noisy data.
   ;
   ; 24-Jun-2015: max(ee) -> median(ee) to prevent a single large error
   ;              bar causing fit to be skipped.
   ;
    IF swtch THEN BEGIN
      j=where(yy NE miss,nj)
     ;
     ; If the number of free parameters is greater than or equal to the
     ; number of data points, then a fit is not possible.
     ;
      IF nj GT n_free THEN BEGIN
       ;
        IF n_tags(template) EQ 0 THEN BEGIN
         ;
         ; TEMPLATE not set - assume single Gaussian fit
         ; ----------------
         ;  - derive initial parameters
         ;
          getmax=max(yy[j],imax)
          init=[max(yy[j])-min(yy[j]),xx[j[imax]],init_wid,min(yy[j]),min(yy[j])]
          IF keyword_set(uniform_backg) THEN init=init[0:3]

         ;
         ; Setting parinfo prevents the fitting routine generating silly fits,
         ; but there is a small risk that interesting events may be missed.
         ;
          parinfo=replicate({fixed: 0, limited: [0,0], limits:[0.d,0.d]},nparams)
         ;
         ; force amplitude to be > 0
          parinfo[0].limited[0]=1
          parinfo[0].limits[0]=0.0
         ; force centroid to be within +/- 0.2 angstroms of initial guess
          parinfo[1].limited=[1,1]
          parinfo[1].limits=[xx[j[imax]]-0.2,xx[j[imax]]+0.2]
         ; force Gaussian width to lie within a fixed range
          parinfo[2].limited=[1,1]
          parinfo[2].limits=wid_range


        ENDIF ELSE BEGIN
         ;
         ; TEMPLATE set 
         ; ------------
         ; If /uniform_backg has been set then init and parinfo need to be
         ; trimmed. 
          init=template.init[0:nparams-1]
          parinfo=parinfo[0:nparams-1]
         ;
         ; Scale the peak values stored in INIT
         ; using the maximum of the current spectrum.
         ;
          i_init=indgen(ngauss)*3
          peaks=init[i_init]
          init[i_init]=peaks/max(peaks)*max(yy[j]) ; changed from yy to yy[j], PRY, 3-May-11
         ;
         ; set guess for background to be minimum of spectrum
         ;
          init[3*ngauss]=min(yy[j])
          IF nback EQ 2 THEN init[3*ngauss+1]=min(yy[j])
        ENDELSE

       ;
       ; This is where the fit is performed
       ; ----------------------------------
        aa = MPFITEXPR(expr, xx[j], yy[j], ee[j], init, $
                       perr=sigmaa, /quiet, bestnorm=bestnorm,yfit=yfit, $
                       parinfo=parinfo,status=status)

        
        IF n_elements(sigmaa) EQ 0 THEN BEGIN
         ;
         ; Sometimes the fit fails to be performed and sigmaa ends up undefined
         ; so I flag this case as bad_pix=5
         ; 
          fitdata.bad_pix[i0,k0]=5b
        ENDIF ELSE BEGIN 

       ;
       ; If parameter tying has been used, then the following computes the
       ; parameter values for the tie parameters.
       ;       
        IF n_tags(template) NE 0 THEN BEGIN
          FOR ig=0,ngauss-1 DO BEGIN
            IF template.lines[ig].peak_tie GE 0 THEN BEGIN
              i_tie=template.lines[ig].peak_tie*3
              aa[ig*3]=aa[i_tie]*template.lines[ig].peak_tie_val
              sigmaa[ig*3]=sigmaa[i_tie]*template.lines[ig].peak_tie_val
            ENDIF
           ;
            IF template.lines[ig].cen_tie GE 0 THEN BEGIN
              i_tie=template.lines[ig].cen_tie*3+1
              aa[ig*3+1]=aa[i_tie]+template.lines[ig].cen_tie_val
              sigmaa[ig*3+1]=sigmaa[i_tie]
            ENDIF
           ;
            IF template.lines[ig].wid_tie GE 0 THEN BEGIN
              i_tie=template.lines[ig].wid_tie*3+2
              aa[ig*3+2]=aa[i_tie]
              sigmaa[ig*3+2]=sigmaa[i_tie]
            ENDIF
          ENDFOR
        ENDIF 

       ;
       ; If any of the errors end up as zero, then I flag the pixel in
       ; 'bad_pixel'. However it's possible that the fits for some of
       ; the lines are perfectly OK so I don't change aa or sigmaa
       ;
       ; PRY, 12-Apr-2011, I've commented out the bad_pix line as I don't see
       ; any use for this information.
       ;
;        chck=where(sigmaa EQ 0.,nchck)
;        IF nchck GT 0 THEN fitdata.bad_pix[i0,k0]=100b+byte(nchck)

          fitdata.aa[*,i0,k0]=aa
          fitdata.sigmaa[*,i0,k0]=sigmaa

          fitdata.x_bg1[i0,k0]=min(xx[j])
          fitdata.x_bg2[i0,k0]=max(xx[j])

          fitdata.chi2[i0,k0]=bestnorm/(n_elements(j) - n_free)

         ;
         ; Compute intensity and intensity error. The formula used for the
         ; intensity error is
         ;   (o_I/I)^2 = 0.5 * ( (o_P/P)^2 + (o_W/W)^2 )
         ; for I=intensity, P=peak and W=width,  o_I, etc. indicates the
         ; 1-sigma error on the quantity. The 0.5 factor assumes that line peak
         ; and width are correlated.
         ;         
          FOR il=0,ngauss-1 DO BEGIN
            i1=il*3    ; line peak
            ic=il*3+1  ; line centroid
            i2=il*3+2  ; line width
            IF keyword_set(perpixel) THEN BEGIN
              getmin=min(abs(wvl-fitdata.aa[ic,i0,k0]),imin)
              IF imin EQ 0 THEN wvl_factor=abs(wvl[1]-wvl[0]) ELSE wvl_factor=abs(wvl[imin]-wvl[imin-1])
            ENDIF ELSE BEGIN
              wvl_factor=1.0
            ENDELSE
           ;
           ; PRY, 13-Feb-2018. In the expression for interr, I had
           ; been dividing by wvl_factor, but int already included
           ; this factor, so the error ended up too large. This is now
           ; fixed.
           ;
            fitdata.int[il,i0,k0]=fitdata.aa[i1,i0,k0]*fitdata.aa[i2,i0,k0]*sqrt(2.*!pi)/wvl_factor
            fitdata.interr[il,i0,k0]=fitdata.int[il,i0,k0]*sqrt(0.5*((fitdata.sigmaa[i1,i0,k0]/fitdata.aa[i1,i0,k0])^2. + (fitdata.sigmaa[i2,i0,k0]/fitdata.aa[i2,i0,k0])^2. ))
          ENDFOR 

        ENDELSE 

        ENDIF ELSE BEGIN
          IF nj EQ 0 THEN fitdata.bad_pix[i0,k0]=3b ELSE fitdata.bad_pix[i0,k0]=2b
        ENDELSE 
      ENDIF ELSE BEGIN
        fitdata.bad_pix[i0,k0]=1b
      ENDELSE 

    IF NOT keyword_set(quiet) THEN BEGIN
      pct=float(ny*(i-xrange[0])+(k-yrange[0])+1)/float(nx)/float(ny)*100.
      progress,pct
    ENDIF

  ENDIF ELSE BEGIN
    fitdata.bad_pix[i0,k0]=4b
  ENDELSE 

  ENDFOR 
  
ENDFOR

IF NOT keyword_set(quiet) THEN BEGIN
  progress,100.,/last
  print,''
ENDIF

;
; The following code sets the REFWVL tag to be the average centroid
; over the raster. Note that missing values are not included in the
; average. 
;
FOR i=0,ngauss-1 DO BEGIN
  cen=reform(fitdata.aa[i*3+1,*,*])
  av_cen=average(cen,missing=0.)
  fitdata.refwvl[i]=av_cen
ENDFOR

;
; Through the fit template it is possible to place limits on the range
; of variability of the fit parameters. E.g., one may want to restrict
; line widths to be between 0.02 and 0.05. The following code checks
; to see if any of these limits were reached, and prints a summary.
;
IF NOT keyword_set(quiet) THEN print,'Checking parameter limits...'
chck=0
FOR i=0,nparams-1 DO BEGIN
  FOR j=0,1 DO BEGIN
    IF parinfo[i].limited[j] EQ 1 THEN BEGIN
      limit=parinfo[i].limits[j]
      arr=reform(fitdata.aa[i,*,*])
      bad_pix=fitdata.bad_pix
      k=where(arr EQ limit AND bad_pix EQ 0,nk)
      IF nk GT 0 THEN BEGIN
        chck=1
        il=i/3
        ip=i-il*3
        CASE ip OF 
          0: pstr='Peak'
          1: pstr='Cent'
          2: pstr='Widt'
        ENDCASE 
        CASE j OF
          0: limstr='Lo'
          1: limstr='Hi'
        ENDCASE 
        IF instrume EQ 'IRIS' THEN wvl_str='f7.2' ELSE wvl_str='f6.2'
        format='(" - Line: ",i2," (",'+wvl_str+',")  -  Parameter: ",a4,"  -  Limit: ",a2,"  -  Pixels: ",i5)'
        print,format=format,il,fitdata.refwvl[il],pstr,limstr,nk
      ENDIF
    ENDIF 
  ENDFOR 
ENDFOR 
IF chck EQ 0 THEN BEGIN
  IF NOT keyword_set(quiet) THEN print,'No limit cases found (this is good!).'
ENDIF ELSE BEGIN
  IF keyword_set(quiet) THEN BEGIN
    print,'%EIS_AUTO_FIT: some parameter limits have been reached. Please run without'
    print,' the /quiet keyword for more details.'
  ENDIF ELSE BEGIN 
    print,'The above parameter limits have been reached. Please consider revising the parameter limits in TEMPLATE.'
  ENDELSE
ENDELSE 


;
; The following prints information about 'bad pixels', i.e., pixels
; stored in the 'bad_pix' array. Fits were not attempted for these
; pixels. 
;
IF NOT keyword_set(quiet) THEN BEGIN
  bad_pix=fitdata.bad_pix
  siz=size(bad_pix)
  n=siz[1]*siz[2]
  chck=where(bad_pix EQ 1,n1)
  chck=where(bad_pix EQ 2,n2)
  chck=where(bad_pix EQ 3,n3)
  chck=where(bad_pix EQ 4,n4)
  IF n1+n2+n3 NE 0 THEN BEGIN
    print,''
    print,'Pixels with no fits: '
    print,'  Too noisy: '+trim(string(format='(f10.1)',float(n1)/float(n)*100.))+'%'
    print,'  Too many missing pixels: '+trim(string(format='(f10.1)',float(n2)/float(n)*100.))+'%'
    print,'  Completely missing spectrum: '+trim(string(format='(f10.1)',float(n3)/float(n)*100.))+'%'
    print,'  Missing offset: '+trim(string(format='(f10.1)',float(n4)/float(n)*100.))+'%'
  ENDIF 
ENDIF 

fitdata.time_stamp=systime()

t2=systime(1)

IF NOT keyword_set(quiet) THEN BEGIN
  print,''
  print,format='("Time taken: ",i5," secs")',round(t2-t1)
  print,''
ENDIF 

END

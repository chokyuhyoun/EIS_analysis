;+
; Function: get_ion_props_young2007
; Description: Returns atomic mass (amu) and Log(T_max) based on 
;              Table 1 of Young et al. (2007), PASJ, 59, S857.
;              https://arxiv.org/pdf/0706.1857
;              Covers Core, Active Region, and Flare lines.
;              Updated values from Chianti v11.0.1
;              https://linelists.chiantidatabase.org/ch_line_list_v11.0.1_150_912.pdf
;-
FUNCTION get_ion_props_young2007, line_id

  ; Normalize input to uppercase for matching
  id = STRUPCASE(STRTRIM(line_id, 2))
  
  mass = 0.
  logt = 0.
  
  ; -------------------------------------------------------
  ; MATCHING LOGIC
  ; Matches substrings common in EIS filenames/headers
  ; -------------------------------------------------------
  
  ; --- HELIUM (Mass ~4) ---
  IF (STRPOS(id, 'HE II') NE -1) THEN BEGIN
      mass = 4.00 & logt = 4.90
  ENDIF
  
  ; --- OXYGEN (Mass ~16) ---
  IF (STRPOS(id, 'O V') NE -1) THEN BEGIN
      mass = 16.00 & logt = 5.35
  ENDIF
  IF (STRPOS(id, 'O VI') NE -1) THEN BEGIN 
      mass = 16.00 & logt = 5.45
  ENDIF

  ; --- MAGNESIUM (Mass ~24.3) ---
  IF (STRPOS(id, 'MG V')   NE -1) THEN BEGIN
      mass = 24.31 & logt = 5.40 
  ENDIF
  IF (STRPOS(id, 'MG VI')  NE -1) THEN BEGIN 
      mass = 24.31 & logt = 5.60 
  ENDIF
  IF (STRPOS(id, 'MG VII') NE -1) THEN BEGIN 
      mass = 24.31 & logt = 5.75 
  ENDIF
  
  ; --- SILICON (Mass ~28.1) ---
  IF (STRPOS(id, 'SI VII') NE -1) THEN BEGIN 
      mass = 28.09 & logt = 5.75 
  ENDIF
  IF (STRPOS(id, 'SI X')   NE -1) THEN BEGIN 
      mass = 28.09 & logt = 6.20 
  ENDIF
  
  ; --- SULFUR (Mass ~32.1) ---
  IF (STRPOS(id, 'S XIII') NE -1) THEN BEGIN 
      mass = 32.06 & logt = 6.45 
  ENDIF

  ; --- NICKEL (Mass ~58.7) ---
  IF (STRPOS(id, 'NI XVII') NE -1) THEN BEGIN
    mass = 58.69 & logt = 6.75
  ENDIF


  ; --- CALCIUM (Mass ~40.1) ---
  IF (STRPOS(id, 'CA XVII') NE -1) THEN BEGIN 
      mass = 40.08 & logt = 6.85 
  ENDIF
  
  ; --- IRON (Mass ~55.85) ---
  IF (STRPOS(id, 'FE VIII') NE -1) THEN BEGIN 
      mass = 55.85 & logt = 5.70 ; Cool loops
  ENDIF
  IF (STRPOS(id, 'FE X')    NE -1) THEN BEGIN 
      mass = 55.85 & logt = 6.05 
  ENDIF
  IF (STRPOS(id, 'FE XI')   NE -1) THEN BEGIN 
      mass = 55.85 & logt = 6.15 
  ENDIF
  IF (STRPOS(id, 'FE XII')  NE -1) THEN BEGIN 
      mass = 55.85 & logt = 6.20 ; Peak EIS sensitivity
  ENDIF
  IF (STRPOS(id, 'FE XIII') NE -1) THEN BEGIN 
      mass = 55.85 & logt = 6.25 
  ENDIF
  IF (STRPOS(id, 'FE XIV')  NE -1) THEN BEGIN 
      mass = 55.85 & logt = 6.35 ; Active Region core
  ENDIF
  IF (STRPOS(id, 'FE XV')   NE -1) THEN BEGIN 
      mass = 55.85 & logt = 6.40 
  ENDIF
  IF (STRPOS(id, 'FE XVI')  NE -1) THEN BEGIN 
      mass = 55.85 & logt = 6.4 
  ENDIF
  IF (STRPOS(id, 'FE XVII') NE -1) THEN BEGIN 
      mass = 55.85 & logt = 6.80 
  ENDIF
  
  ; --- FLARE LINES ---
  IF (STRPOS(id, 'FE XXIII') NE -1) THEN BEGIN 
      mass = 55.85 & logt = 7.15 
  ENDIF
  IF (STRPOS(id, 'FE XXIV')  NE -1) THEN BEGIN 
      mass = 55.85 & logt = 7.20 
  ENDIF

  RETURN, {line_id:line_id, mass: mass, logt: logt}
END

;+
; Main Procedure
;-
PRO analyze_eis, path

  f = file_search(path+'/eis_l1_*.fits', count = n)
  if n eq 0 then begin 
    cd, path
    f = file_search(path+'/eis_l0_*.fits', count = n)
    eis_prep, f, /def, /save
    f = file_search(path+'/eis_l1_*.fits', count = n)
    cd, '..'
  endif 
;  stop
  print, 'Starting ', f
  ; Constants
  k_boltz   = 1.38e-16  ; erg/K
  c_light   = 2.9979e5  ; km/s
  amu       = 1.66e-24  ; g
  w_inst    = 0.056     ; Instrumental Width (Angstroms, approx for 1" slit)
  fwhm_to_1e  = 1.0 / (2.0 * SQRT(ALOG(2.0))) 
      
  ; 3. ITERATE WINDOWS
  hdr = fitshead2struct(HEADFITS(f))
  n_win = hdr.nwin
  PRINT, 'Processing ', n_win, ' spectral windows...'
  int_map_all = !null
  vel_map_all = !null
  v_nt_map_all = !null
  ion_data_list = !null
  fitdata_list = list()
  FOR i = 0, n_win - 1 DO BEGIN
stop
      ; --- A. Load Data ---
      ; /quiet suppresses excessive text output
      wd = eis_getwindata(f, i, /refill, /quiet)
      IF (SIZE(wd, /TYPE) NE 8) THEN CONTINUE

      line_id = wd.line_id
      
      ; --- B. Identify Ion & Params ---
      ; Use the helper function based on Young et al. 2007
      ion_data = get_ion_props_young2007(line_id)
      ion_data = create_struct(ion_data, 'win', i)
      temp_k   = 10.^(ion_data.logt)
      mass_g   = ion_data.mass * amu

      PRINT, FORMAT='("Window: ",I02," | Line: ",A-15," | LogT: ",F4.2," | Mass: ",F5.2)', $
             i, line_id, ion_data.logt, ion_data.mass
      mm = 1
      ; --- C. Gaussian Fit ---
      eis_auto_fit_, wd, fitdata, /quiet ; original function "eis_auto_fit" has problem with "nfree" variable. I changed 
      eis_fit_viewer, wd, fitdata
      ;fitdata.refwvl : ch_line_list_v11.0.1_150_912.pdf
      ;ion_data
      ;eis_wvl_select, wd, wvl_select
      ;eis_auto_fit_, wd, fitdata, wvl_select=wvl_select, /quiet

stop
      ; --- D. Extract Maps ---
      if mm then begin
        int_map = eis_get_fitdata(fitdata, /int, /map)
        int_map.data[where(int_map.data eq int_map.missing, /null)] = !values.f_nan
        int_map.missing = !values.f_nan
        
        vel_map = eis_get_fitdata(fitdata, /vel, /map)
        vel_map.data[where(vel_map.data eq vel_map.missing, /null)] = !values.f_nan
        vel_map.missing = !values.f_nan
        
        v_nt_map = eis_get_fitdata(fitdata, /wid, /map, THERMAL_WID=[ion_data.logt, ion_data.mass]) ; FWHM in Angstroms, instr. & thermal width subtracted
        v_nt_map.data[where(v_nt_map.data eq v_nt_map.missing, /null)] = !values.f_nan 
  
        ; 2. Convert Source FWHM (Angstroms) to 1/e Doppler Width (km/s)
        v_nt_map.data  = (v_nt_map.data / fitdata.refwvl[0]) * c_light * fwhm_to_1e
        v_nt_map.missing = !values.f_nan
  
        int_map_all = [int_map_all, int_map]
        vel_map_all = [vel_map_all, vel_map]
        v_nt_map_all = [v_nt_map_all, v_nt_map]
        ion_data_list = [ion_data_list, ion_data]
        fitdata_list.add, fitdata
      endif
stop
  ENDFOR
  save, int_map_all, vel_map_all, v_nt_map_all, ion_data_list, fitdata_list, $
    filename=path+path_sep()+'EIS_fit_results.sav'
END

path = '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff'
cd, path
path_list = file_search(path, '*', /test_dir, /full)
;analyze_eis, path_list[14]
;analyze_eis, '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/20190412_145240'
analyze_eis, '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/target1'
end

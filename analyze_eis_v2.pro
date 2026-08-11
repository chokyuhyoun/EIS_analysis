all_ions  = ['Fe XI', 'Fe XI', 'Fe X', 'Fe VIII', 'Fe XII', $
             'Fe XII', 'Fe XII', 'Fe XII', 'Fe XIII', 'Fe XIII', $
             'Fe XIII', 'Fe XIV', 'Ni XVII', 'Fe XIII', 'He II', $
             'Si XIII', 'Si X', 'Si X', 'Fe XVI', 'S X', $
             'Fe XIV', 'Fe XIV', 'Fe XIV', 'Si VII', 'Mg VII', $
             'Fe XV', 'Ca XVII']
all_waves = [180.4010d, 182.1670d, 184.5370d, 185.2130d, 186.8870d, $
             192.3940d, 193.5090d, 195.1190d, 200.0210d, 201.1260d, $
             202.0440d, 211.3170d, 249.1890d, 251.9520d, 256.3180d, $
             256.6850d, 258.3740d, 261.0560d, 262.9760d, 264.2300d, $
             264.7880d, 270.5200d, 274.2030d, 275.3680d, 278.3940d, $
             284.1630d, 192.8200d]
all_logts = [6.15, 6.15, 6.05, 5.70, 6.20, $
             6.20, 6.20, 6.20, 6.25, 6.25, $
             6.30, 6.35, 6.75, 6.25, 4.90, $
             6.45, 6.15, 6.20, 6.80, 6.25, $
             6.30, 6.35, 6.35, 5.75, 5.75, $
             6.40, 6.70]
; Ca XVII 192.820 / 6.7 --> blended with O V 192.9040 / 5.35
all_masses = [55.845, 55.845, 55.845, 55.845, 55.845, 55.845, 55.845, 55.845, 55.845, 55.845, $
              55.845, 55.845, 58.693, 55.845, 4.0026, 28.085, 28.085, 28.085, 55.845, 32.060, $
              55.845, 55.845, 55.845, 28.085, 24.305, 55.845, 40.078]

n_elements = n_elements(all_ions)
line_data0 = replicate({win: 0, ion: '', wave: 0.0d, logt: 0.0, mass: 0.0}, n_elements)

for i = 0, n_elements - 1 do begin &$
    line_data0[i].ion  = all_ions[i] &$
    line_data0[i].wave = all_waves[i] &$
    line_data0[i].logt = all_logts[i] &$
    line_data0[i].mass = all_masses[i] &$
endfor


path = '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/3791'
cd, path
f = file_search(path+'/eis_l1_*.fits', count = n)
if n eq 0 then begin
  f = file_search(path+'/eis_l0_*.fits', count = n)
  eis_prep, f, /def, /save
  f = file_search(path+'/eis_l1_*.fits', count = n)
endif


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
line_data = !null
fitdata_list = list()


yr_for_wvcor = [450, 511] ; in pix

FOR i = 0, n_win - 1 DO BEGIN
  wd = eis_getwindata(f, i, /refill, /quiet)
  loc = where((line_data0.wave gt min(wd.wvl)) and $
              (line_data0.wave lt max(wd.wvl)), /null)
  for j=0, n_elements(loc)-1 do begin
    dum = line_data0[loc[j]]
    dum.win = i
    PRINT, FORMAT='("Window: ",I02," | Wave: ", F8.3, " | Line: ",A-6," | LogT: ",F4.2," | Mass: ",F5.2)', $
      dum.win, dum.wave, dum.ion, dum.logt, dum.mass
    eis_auto_fit_, wd, fitdata, /quiet, wvlpix=dum.wave + 0.2*[-1, 1], /uniform_backg
    fitdata.refwvl =  dum.wave
    if n_elements(yr_for_wvcor) ne 0 then begin
      newfitdata = eis_update_fitdata(fitdata, yrange=yr_for_wvcor, offset=offset)
      eis_auto_fit_, wd, fitdata, /quiet, wvlpix=dum.wave + 0.2*[-1, 1], /uniform_backg, $
                     offset=offset
    endif
;    eis_fit_viewer, wd, fitdata

    int_map = eis_get_fitdata(fitdata, /int, /map)
    int_map.data[where(int_map.data eq int_map.missing, /null)] = !values.f_nan
    int_map.missing = !values.f_nan
    int_map = create_struct(int_map, dum)

    vel_map = eis_get_fitdata(fitdata, /vel, /map)
    vel_map.data[where(vel_map.data eq vel_map.missing, /null)] = !values.f_nan
    vel_map.missing = !values.f_nan
    vel_map = create_struct(vel_map, dum)

    v_nt_map = eis_get_fitdata(fitdata, /wid, /map, THERMAL_WID=[dum.logt, dum.mass]) ; FWHM in Angstroms, instr. & thermal width subtracted
    v_nt_map.data[where(v_nt_map.data eq v_nt_map.missing, /null)] = !values.f_nan

    ; 2. Convert Source FWHM (Angstroms) to 1/e Doppler Width (km/s)
    v_nt_map.data  = (v_nt_map.data / fitdata.refwvl[0]) * c_light * fwhm_to_1e
    v_nt_map.missing = !values.f_nan
    v_nt_map = create_struct(v_nt_map, dum)

    int_map_all = [int_map_all, int_map]
    vel_map_all = [vel_map_all, vel_map]
    v_nt_map_all = [v_nt_map_all, v_nt_map]
    line_data = [line_data, dum]
    fitdata_list.add, fitdata
  endfor
endfor
save_name = (n_elements(yr_for_wvcor) eq 0) ? 'EIS_fit_results.sav' : 'EIS_fit_wvcor_results.sav'
save, int_map_all, vel_map_all, v_nt_map_all, line_data, fitdata_list, $
  filename=path+path_sep()+save_name

end    
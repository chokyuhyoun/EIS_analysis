; align between Hinode and AIA
; 
;
;path = '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/20190412_145240'
;restore, path+'/instr_info.sav'
;if 0 then begin
;  iris_path = (((iris_info.umodes).split('/'))[0:7]).join('/')
;  iris_file = file_search(iris_path+'/*raster*.fits')
;  get_iris_si_iv_1403_fit, iris_file, amp, vel, nth, iris_xp, iris_yp
;  save, iris_file, amp, vel, nth, iris_xp, iris_yp, $
;    filename=path+'/iris_fit_res.sav'
;endif else begin
;  restore, path+'/iris_fit_res.sav'
;endelse
;
restore, path+'/EIS_fit_results.sav', /verbose
; int_map_all, vel_map_all, v_nt_map_all, ion_data_list

eis_mid_time = 0.5*(anytim(eis_info.date_obs) + anytim(eis_info.date_end))
loc195 = (where(strmatch(ion_data_list.line_id, '*195*')))[0]
wave195 = float(strmid(ion_data_list[loc195].line_id, $
                       strlen(ion_data_list[loc195].line_id)-7))
offset_y = eis_ccd_offset(wave195)
get_xp_yp, int_map_all[loc195], eis_xp0, eis_yp0
eis_yp0 -= offset_y[0]
eis195 = int_map_all[loc195].data
;
;eis_file = file_search(path+'/*eis*l1*.fits')
;wd = eis_getwindata(eis_file[0], loc195, /refill, /quiet)

iris_mid_time = anytim(0.5*(anytim(iris_info.starttime)+anytim(iris_info.stoptime)), /ccsds)
dum = min(abs(eis_xp - iris_info.xcen), eis_close_step)
eis_time_diff = wd.time - wd.time[eis_close_step]

;ssw_jsoc_time2data, wd.time_ccsds[eis_close_step], wd.time_ccsds[eis_close_step-1], $
;  aia193_index, ds='aia.lev1_euv_12s', wave=193, $
;  /files_only, aia193_file, /silent
;aia193_file = aia193_file[0]
;aia193_index = aia193_index[0]
;read_sdo, aia193_file, dum, aia193_data, /noshell, /use_shared_lib
;get_xp_yp, aia193_index, aia193_xp, aia193_yp

eis_xp = eis_xp0
eis_yp = eis_yp0
spix = [0., 0]

for j=0, 20 do begin
  eis_xp = eis_xp0 + spix[0]*int_map_all[0].dx
  eis_yp = eis_yp0 + spix[1]*int_map_all[0].dy
  eis_nxp0 = !null
  for i=0, n_elements(eis_xp)-1 do $
    eis_nxp0 = [eis_nxp0, (rot_xy(eis_xp[i], mean(eis_yp), -eis_time_diff[i]))[0]]
  
  aia193_eis_npixx = interpol(findgen(n_elements(aia193_xp)), aia193_xp, eis_nxp0)
  aia193_eis_pixy = interpol(findgen(n_elements(aia193_yp)), aia193_yp, eis_yp)
  aia193_npart = interpolate(aia193_data, aia193_eis_npixx+spix[0], aia193_eis_pixy+spix[1], /grid)
  spix = alignoffset(aia193_npart, eis195) + spix
  print, spix
endfor
eis_nxp0 = !null
for i=0, n_elements(eis_xp)-1 do $
  eis_nxp0 = [eis_nxp0, (rot_xy(eis_xp[i], mean(eis_yp), -eis_time_diff[i]))[0]]
eis_nxp = eis_nxp0 + spix[0]*int_map_all[0].dx
eis_nyp = eis_yp + spix[1]*int_map_all[0].dy
end
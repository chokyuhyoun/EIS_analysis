nlfff_path = '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff'
f = file_search(nlfff_path, 'nlfff_results.sav')

;for i=0, n_elements(f)-1 do begin
i = 14
base_path = file_dirname(f[i])

cd, base_path
print, f[i]
restore, f[i]
restore, base_path+'/instr_info.sav', /verbose
; iris_info, eis_info, index_blos

t0 = index_i0.date_obs
t1 = anytim(reltime(index_i0.date_obs, min=10), /ccsds)
fov = 50

ssw_jsoc_time2data, t0, t1, $
  index_blos, ds='hmi.M_45s', segment='magnetogram', $
  files_only=1, filenames, /silent
index_blos = index_blos[0]
read_sdo, filenames[0], dum, data_blos, /sil, /noshell, /use_shared_lib
get_xp_yp, index_blos, xp_blos, yp_blos, data=data_blos

ssw_jsoc_time2data, t0, t1, $
  index_aia0, ds='aia.lev1_euv_12s', wave=304, $
  files_only=1, filenames, /silent
index_aia0 = index_aia0[0]
read_sdo, filenames[0], dum, data_aia0, /sil, /noshell, /use_shared_lib
get_xp_yp, index_aia0, xp_aia0, yp_aia0

restore, base_path+'/EIS_fit_results.sav'
; int_map_all, vel_map_all, v_nt_map_all, ion_data_list
ind = where(strmatch(ion_data_list.line_id, 'He II 256*'))
int_map = int_map_all[ind]

aia_lct, rr, gg, bb, wave=304
w01 = window(dim=[8e2, 8e2])
im01 = image_(data_aia0, xp_aia0, yp_aia0, /current, $
              min=0, max=1e2, rgb_table=[[rr], [gg], [bb]], $
              xr=int_map.xc+[-1, 1]*fov, yr=int_map.yc+[-1, 1]*fov)
im02 = image_(int_map, over=im01)
im03 = image_(data_blos, xp_blos, yp_blos, over=im01, min=-500, max=500)
ng_blink, [im01, im02, im03]

;====================================

restore, base_path+'/iris_fit_res.sav'
;amp, vel, nth, iris_xp, iris_yp

st = iris_info.starttime
et = iris_info.stoptime
mt0 = anytim(0.5*(anytim(st) + anytim(et)), /ccsds)
mt1 = anytim(0.5*(anytim(st) + anytim(et)) + 24, /ccsds)
;ssw_jsoc_time2data, mt0, mt1, $
ssw_jsoc_time2data, t0, t1, $  
  index_aia1, ds='aia.lev1_uv_24s', wave=1600, $
  files_only=1, filenames, /silent
index_aia1 = index_aia1[0]
read_sdo, filenames[0], dum, data_aia1, /sil, /noshell, /use_shared_lib
get_xp_yp, index_aia1, xp_aia1, yp_aia1

aia_lct, rr, gg, bb, wave=1600
w02 = window(dim=[8e2, 8e2])
im11 = image_(data_aia1, xp_aia1, yp_aia1, /current, $
              min=0, max=5e2, rgb_table=[[rr], [gg], [bb]], $
              xr=int_map.xc+[-1, 1]*fov, yr=int_map.yc+[-1, 1]*fov)
im12 = image_(amp, iris_xp, iris_yp, over=im11, min=0, max=5e2)
im13 = image_(data_blos, xp_blos, yp_blos, over=im11, min=-500, max=500)
ng_blink, [im11, im12, im13]
end
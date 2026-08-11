FUNCTION sync_fov_wheel, oWin, x, y, delta, keymods

  graphic_obj = oWin.HitTest(x, y)
  IF ~ISA(graphic_obj) THEN RETURN, -1
  if n_elements(graphic_obj) ne 1 then return, -1
  xy = graphic_obj.ConvertCoord(x, y, /DEVICE, /TO_DATA)
  oWin.refresh, /disable
  u = oWin.UVALUE
  IF delta GT 0 THEN zoom = 0.90D ELSE zoom = 1.1D
  u.xr = (u.xr-xy[0])*zoom + xy[0]
  u.yr = (u.yr-xy[1])*zoom + xy[1]
  for i=0, n_elements(u.ims)-1 do begin
    u.ims[i].xr = u.xr
    u.ims[i].yr = u.yr
  endfor
  oWin.uvalue = u
  oWin.Refresh
  RETURN, 0
END

function mousedown, oWin, x, y, iButton, KeyMods, nClicks
  graphic_obj = oWin.HitTest(x, y)
  IF ~ISA(graphic_obj) THEN RETURN, -1

  u = oWin.uvalue
  if ibutton eq 2 then begin
    u.midbuttonDown = 1
    u.midcurx = x
    u.midcury = y
  endif
  if ibutton eq 1 then begin
    oWin.refresh, /disable
    xy_data = (oWin.ConvertCoord(x, y, /DEVICE, /TO_DATA))[0:1]
    if ~finite(u.leftcurx[0]) then begin
      u.leftcurx = list(xy_data[0])
      u.leftcury = list(xy_data[1])
    endif else begin
      u.leftcurx.add, xy_data[0]
      u.leftcury.add, xy_data[1]
    endelse
    for k=0, n_elements(u.point_plot)-1 do begin
      u.point_plot[k].setdata, u.leftcurx.toarray(), u.leftcury.toarray()
    endfor
    oWin.refresh
  endif
  oWin.uvalue = u
  return, 0
end

function del_points, oWin
  u = oWin.uvalue
  u.leftcurx = list(!values.f_nan)
  u.leftcury = list(!values.f_nan)
  for k=0, n_elements(u.point_plot)-1 do begin
    u.point_plot[k].setdata, u.leftcurx.toarray(), u.leftcury.toarray()
  endfor    
  oWin.uvalue = u
  return, 0
end

function hide_points, oWin
  u = oWin.uvalue
  for k=0, n_elements(oWin.uvalue.point_plot)-1 do begin
    u.point_plot[k].hide = (u.point_plot[k].hide + 1) mod 2 
  endfor
  oWin.uvalue = u
  return, 0
end  
 
function sync_fov_move, oWin, x, y, button
  if oWin.uvalue.midbuttondown eq 0 then return, 0
  graphic_obj = oWin.HitTest(x, y)
  IF ~ISA(graphic_obj) THEN RETURN, -1
  if n_elements(graphic_obj) ne 1 then return, -1
  oWin.refresh, /disable
  u = oWin.uvalue
  xy = graphic_obj.ConvertCoord([u.midcurx, x], [u.midcury, y], /DEVICE, /TO_DATA)
  u.xr -= xy[0, 1] - xy[0, 0]
  u.yr -= xy[1, 1] - xy[1, 0]
  for i=0, n_elements(u.ims)-1 do begin
    u.ims[i].xr = u.xr
    u.ims[i].yr = u.yr
  endfor
  u.midcurx = x
  u.midcury = y
  oWin.uvalue = u
  oWin.Refresh
  RETURN, 0
END

function mouseup, oWin, x, y, iButton, KeyMods, nClicks
  u = oWin.uvalue
  if ibutton eq 2 then u.midbuttondown = 0
  oWin.uvalue = u
  return, 0
end

function instr_rotation, window, i
  window.refresh, /disable
  u = window.uvalue
  u.i = i
  for k=0, n_elements(u.ims)-1 do u.ims[k].hide = 1
  for k=0, 2 do begin
    u.ims[k, i].hide = 0
    u.ims[k, i].order, /bring_to_front
  endfor
  for k=0, n_elements(u.point_plot)-1 do begin
    u.point_plot[k].order, /bring_to_front
  endfor
  u.ims[0, 0].title = u.im01_titles[i]
  window.uvalue = u
  window.refresh 
  return, 0
end
 
function eis_rotation, window, j
  window.refresh, /disable
  u = window.uvalue
  u.j = j
  int_min = u.ims[0, 0].min
  vel_min = u.ims[1, 0].min
  vel_max = u.ims[1, 0].max
  nth_min = u.ims[2, 0].min
  nth_max = u.ims[2, 0].max
  eis_yp = u.eis_yp0 - u.offset_y[j]
  u.ims[0, 0].setdata, u.int_map[j].data, u.eis_xp, eis_yp
  u.ims[1, 0].setdata, u.vel_map[j].data, u.eis_xp, eis_yp
  u.ims[2, 0].setdata, u.v_nt_map[j].data, u.eis_xp, eis_yp
  hist = histogram(u.int_map[j].data, loc=xbin)
  cdf = total(hist, /cum)/total(finite(u.int_map[j].data))
  int_max = xbin[min(where(cdf gt 0.98))]
  u.ims[0, 0].min = 0
  u.ims[0, 0].max = int_max
  u.ims[1, 0].min = vel_min
  u.ims[1, 0].max = vel_max
  u.ims[2, 0].min = nth_min
  u.ims[2, 0].max = nth_max
  u.im01_titles[0] = 'Int: ' + u.int_map[u.j].ion + ' ' + $
    string(u.int_map[u.j].wave, f='(f7.3)') + ' $\AA$' + $
    ' (log T='+string(u.int_map[u.j].logt, f='(f0.1)')+' K)'
  if u.i eq 0 then u.ims[0, 0].title = u.im01_titles[u.i]
  window.uvalue = u
  window.refresh   
  return, 0
end

function blink, window, $
  IsASCII, Character, KeyValue, X, Y, Press, Release, KeyMods
  window.refresh, /disable
  IF release THEN RETURN, 1
  char = string(character)
  if isASCII eq 0 then begin
    if (keyvalue eq 6) or (keyvalue eq 5) then begin
      del = (keyvalue eq 6) ? 1 : -1
      j = (window.uvalue.j + del + n_elements(window.uvalue.int_map)) $
            mod n_elements(window.uvalue.int_map)
      dum = eis_rotation(window, j)
    endif      
  endif else begin
    if character ge 49 and character le 57 then begin
      num_str = string(indgen((size(window.uvalue.ims))[2])+1, f='(i0)')
      if total(strmatch(num_str, char) eq 1) then begin
        i = float(char)
        dum = instr_rotation(window, i-1)
      endif
    endif
    if char eq 'q' then begin 
      window.close
      return, 0
    endif
    if char eq 'd' then dum = del_points(window)
    if char eq 'h' then dum = hide_points(window)
  endelse
  window.refresh 
  return, 0
end

;; analysis
;

path = '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/3791'
fig_save = 0
restore, path+'/instr_info.sav';, /verbose
if n_elements(iris_info) gt 1 then begin ; among the IRIS obs, select largest xfov
  xfov = !null
  for i=0, n_elements(iris_info)-1 do xfov = [xfov, iris_info[i].xfov]  
  iris_info = (iris_info[where(xfov eq max(xfov))])[0]
endif

restore, path+'/EIS_fit_wvcor_results.sav';, /verbose
;restore, path+'/EIS_fit_results.sav';, /verbose
; int_map_all, vel_map_all, v_nt_map_all, line_data, fitdata_list

f = file_search(path+'/si_iv_fit_res.sav', count=n)
if n eq 0 then begin
  iris_path = (((iris_info.umodes).split('/'))[0:7]).join('/')
  iris_file = file_search(iris_path+'/*raster*.fits')

;  iris_file = iris_file[2]
  if n_elements(iris_file) gt 1 then begin &$
    hdr = !null  &$
    for k=0, n_elements(iris_file)-1 do begin &$
      hdr = [hdr, fitshead2struct(headfits(iris_file[k]))] &$
    endfor &$
    dum = min(abs(anytim(hdr.date_obs) - anytim(int_map_all[0].time)), ind) &$ 
    iris_file = iris_file[ind] &$
  endif
  si_iv_fit_res = iris_si_iv_fit(iris_file)
  save, si_iv_fit_res, filename=path+'/si_iv_fit_res.sav'
endif else begin
  restore, path+'/si_iv_fit_res.sav'
endelse

read_sdo, blos_filename, dum, data_blos, /sil, /noshell, /use_shared_lib
get_xp_yp, index_blos[0], hmi_xp, hmi_yp, hmi_xxp, hmi_yyp, data=data_blos

;restore, path+'/nlfff_results.sav'
;iris_xpix = interpol(findgen(n_elements(res_xp)), res_xp, iris_xp)
;iris_ypix = interpol(findgen(n_elements(res_yp)), res_yp, iris_yp)
;;stop
;cos_i_iris = interpolate(reform(abs(cos_i0[*, *, 6])), iris_xpix, iris_ypix, $
;                         missing=!values.f_nan)

dum = where(int_map_all.mass ne 0, /null)
order = sort(int_map_all[dum].logt)
int_map = (int_map_all[dum])[order]
vel_map = (vel_map_all[dum])[order]
v_nt_map = (v_nt_map_all[dum])[order]
get_xp_yp, int_map[0], eis_xp, eis_yp0

i = 0 ; instr
j = 0 ; EIS line id
loc = 0
w03 = window(dim=[1e3, 5e2])
w03.refresh, /disable
hist = histogram(int_map[j].data, loc=xbin)
pdf = total(hist, /cum)/total(finite(int_map[j].data))
int_max = xbin[min(where(pdf gt 0.99))]
xmar = 0.07
xs = (1.-6*xmar)/3.
xi = xmar + findgen(3)*(xs+2*xmar)
yi = 0.15
yf = 0.9

offset_y = eis_ccd_offset(int_map.wave)
eis_yp = eis_yp0 - offset_y[i]
xr = minmax(si_iv_fit_res.xp) + 30.*[-1, 1]
yr = minmax(si_iv_fit_res.yp) + 30.*[-1, 1]
im01 = IMAGE_(int_map[j].data, eis_xp, eis_yp, /current, aspect_r=1, $
              pos=[xi[0], yi, xi[0]+xs, yf], $
              min=0, max=int_max, cb=cb01, xtickdir=1, ytickdir=1, $
              xr=xr, yr=yr, xtickinterval=50, $
              xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)')
im02 = IMAGE_(vel_map[j].data, eis_xp, eis_yp, /current, aspect_r=1, $
              pos=[xi[1], yi, xi[1]+xs, yf], $
              TITLE='Vel', $
              rgb_table=33, min=-30, max=30, cb=cb02, xtickdir=1, ytickdir=1, $
              xr=xr, yr=yr, xtickinterval=50, $
              xtitle='Solar X (arcsec)')
cb02.title = 'km/s'
im03 = IMAGE_(v_nt_map[j].data, eis_xp, eis_yp, /current, aspect_r=1, $
              pos=[xi[2], yi, xi[2]+xs, yf], $
              TITLE='V_nt', $
              rgb_table=4, min=0, max=60, cb=cb03, xtickdir=1, ytickdir=1, $
              xr=xr, yr=yr, xtickinterval=50, $
              xtitle='Solar X (arcsec)')
cb03.title = 'km/s'

im011 = image_(si_iv_fit_res.amp, si_iv_fit_res.xp, si_iv_fit_res.yp, $
               over=im01, min=0, max=4e2)
im021 = image_(si_iv_fit_res.vel, si_iv_fit_res.xp, si_iv_fit_res.yp, $
               over=im02, min=-3e1, max=3e1, rgb_table=33)
im031 = image_(si_iv_fit_res.nth, si_iv_fit_res.xp, si_iv_fit_res.yp, $
               over=im03, min=0, max=3e1, rgb_table=4)
;height_pix = 5
;im013 = image_(cos_i0[*, *, height_pix], res_xp, res_yp, over=im01, min=0, max=1, rgb_table=22)
;im023 = image_(cos_i0[*, *, height_pix], res_xp, res_yp, over=im02, min=0, max=1, rgb_table=22)
;im033 = image_(cos_i0[*, *, height_pix], res_xp, res_yp, over=im03, min=0, max=1, rgb_table=22)
;im013.title = 'Cos i (B vs. LOS)'
;nlfff_ims = [im013, im023, im033]

im012 = image_(data_blos, hmi_xp, hmi_yp, over=im01, min=-500, max=500, rgb_table=0)
im022 = image_(data_blos, hmi_xp, hmi_yp, over=im02, min=-500, max=500, rgb_table=0)
im032 = image_(data_blos, hmi_xp, hmi_yp, over=im03, min=-500, max=500, rgb_table=0)

ssw_jsoc_time2data, int_map_all[0].time, reltime(int_map_all[0].time, min=0.2), $
  index_aia1, ds='aia.lev1_euv_12s', wave=171, $
  files_only=1, filenames, /silent
read_sdo, filenames[0], dum, data_aia, /sil, /use_shared_lib
get_xp_yp, index_aia1[0], aia_xp, aia_yp

aia_lct, rr, gg, bb, wave=index_aia1[0].wavelnth
ct_aia = [[rr], [gg], [bb]]
im013 = image_(data_aia, aia_xp, aia_yp, over=im01, min=0, max=3e3, rgb_table=ct_aia)
im023 = image_(data_aia, aia_xp, aia_yp, over=im02, min=0, max=3e3, rgb_table=ct_aia)
im033 = image_(data_aia, aia_xp, aia_yp, over=im03, min=0, max=3e3, rgb_table=ct_aia)

ssw_jsoc_time2data, int_map_all[0].time, reltime(int_map_all[0].time, min=0.4), $
  index_aia1, ds='aia.lev1_euv_12s', wave=304, $
  files_only=1, filenames, /silent
read_sdo, filenames[0], dum, data_aia, /sil, /use_shared_lib
get_xp_yp, index_aia1[0], aia_xp, aia_yp

aia_lct, rr, gg, bb, wave=index_aia1[0].wavelnth
ct_aia = [[rr], [gg], [bb]]
im014 = image_(alog10(data_aia), aia_xp, aia_yp, over=im01, min=0, max=2, rgb_table=ct_aia)
im024 = image_(alog10(data_aia), aia_xp, aia_yp, over=im02, min=0, max=2, rgb_table=ct_aia)
im034 = image_(alog10(data_aia), aia_xp, aia_yp, over=im03, min=0, max=2, rgb_table=ct_aia)

ssw_jsoc_time2data, int_map_all[0].time, reltime(int_map_all[0].time, min=0.2), $
  index_aia1, ds='aia.lev1_euv_12s', wave=193, $
  files_only=1, filenames, /silent
read_sdo, filenames[0], dum, data_aia, /sil, /use_shared_lib
get_xp_yp, index_aia1[0], aia_xp, aia_yp

aia_lct, rr, gg, bb, wave=index_aia1[0].wavelnth
ct_aia = [[rr], [gg], [bb]]
im015 = image_(data_aia, aia_xp, aia_yp, over=im01, min=0, max=3e3, rgb_table=ct_aia)
im025 = image_(data_aia, aia_xp, aia_yp, over=im02, min=0, max=3e3, rgb_table=ct_aia)
im035 = image_(data_aia, aia_xp, aia_yp, over=im03, min=0, max=3e3, rgb_table=ct_aia)

t01 = text(0.03, 0.045, 'EIS  obs time: '+eis_str.date_obs+' -- '+eis_str.date_end)
t02 = text(0.03, 0.02, 'IRIS obs time: '+si_iv_fit_res.header.date_obs+' -- '+$
                                         si_iv_fit_res.header.date_end)

ims = [[im01, im02, im03], $ ; EIS
       [im011, im021, im031], $ ; IRIS
       [im012, im022, im032], $ ; b_los
       [im013, im023, im033], $ ; AIA 171  
       [im014, im024, im034], $ ; AIA 304
       [im015, im025, im035]]   ; AIA 193
im01_titles = ['Int: ' + int_map[j].ion + ' ' + $
               string(int_map[j].wave, f='(f7.3)') + ' $\AA$' + $
               ' (log T='+string(int_map[j].logt, f='(f0.1)')+' K)', $ 
               'Int: '+si_iv_fit_res.line_id+' $\AA$ (log T=4.8 K)', $
               'B$_{LOS}$', $
               'AIA 171 $\AA$', $ 
               'AIA 304 $\AA$', $
               'AIA 193 $\AA$']
leftcurx = list(!values.f_nan)
leftcury = list(!values.f_nan)
p01 = plot(leftcurx.toarray(), leftcury.toarray(), '+', over=im01)                
p02 = plot(leftcurx.toarray(), leftcury.toarray(), '+', over=im02)
p03 = plot(leftcurx.toarray(), leftcury.toarray(), '+', over=im03)
p04 = plot(leftcurx.toarray(), leftcury.toarray(), 'Xw', over=im01)
p05 = plot(leftcurx.toarray(), leftcury.toarray(), 'Xw', over=im02)
p06 = plot(leftcurx.toarray(), leftcury.toarray(), 'Xw', over=im03)
point_plot = [p01, p02, p03, p04, p05, p06]

w03.uvalue = {ims:ims, im01_titles:im01_titles, $
;              nlfff_ims: [im013, im023, im033], $
              eis_xp:eis_xp, eis_yp0:eis_yp0, $
              int_map:int_map, vel_map:vel_map, v_nt_map:v_nt_map, $
;              cos_i_iris:cos_i_iris, $
;              iris_xpix:iris_xpix, iris_ypix:iris_ypix, $
              si_iv_fit_res:si_iv_fit_res, i:i, j:j, $
              loc:loc, offset_y:offset_y, xr:xr, yr:yr, $
              midbuttondown:0L, midcurx:0l, midcury:0l, $
              leftcurx:leftcurx, leftcury:leftcury, $
              point_plot:point_plot}
w03.keyboard_handler='blink'
w03.MOUSE_WHEEL_HANDLER = 'sync_fov_wheel'
w03.mouse_down_handler = 'mousedown'
w03.mouse_motion_handler = 'sync_fov_move'
w03.mouse_up_handler = 'mouseup'
w03.refresh

if fig_save then begin
  file_mkdir, path+'/figures'
  dum = instr_rotation(w03, 2)
  w03.save, path+'/figures/00 HMI_B_los.png', resol=200
  dum = instr_rotation(w03, 1)
  w03.save, path+'/figures/01 IRIS '+si_iv_fit_res.line_id+'.png', resol=200
  dum = instr_rotation(w03, 3)
  w03.save, path+'/figures/02 AIA 171.png', resol=200
  dum = instr_rotation(w03, 0)
  for k=0, n_elements(int_map)-1 do begin &$
    dum = eis_rotation(w03, k) &$
    filename = string(format='(a, "/figures/",i02," EIS ", a, " ", f7.3, ".png")' , $
                      path, k+3, int_map[k].ion, int_map[k].wave) &$
    w03.save, filename, resol=200 &$
  endfor
  dum = eis_rotation(w03, 0)
endif
end
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
  objs = [u.ims, u.iris_ims, u.nlfff_ims]
  for i=0, n_elements(objs)-1 do begin
    objs[i].xr = u.xr
    objs[i].yr = u.yr
  endfor
  oWin.uvalue = u
  oWin.Refresh
  RETURN, 0
END

function mousedown, oWin, x, y, iButton, KeyMods, nClicks
  if ibutton ne 2 then return, 0
  u = oWin.uvalue
  u.buttonDown = 1
  u.curx = x
  u.cury = y
  oWin.uvalue = u
  return, 0
end
  
function sync_fov_move, oWin, x, y, button
  if oWin.uvalue.buttondown eq 0 then return, 0
  graphic_obj = oWin.HitTest(x, y)
  IF ~ISA(graphic_obj) THEN RETURN, -1
  if n_elements(graphic_obj) ne 1 then return, -1
  oWin.refresh, /disable
  u = oWin.uvalue
  xy = graphic_obj.ConvertCoord([u.curx, x], [u.cury, y], /DEVICE, /TO_DATA)
  u.xr -= xy[0, 1] - xy[0, 0]
  u.yr -= xy[1, 1] - xy[1, 0]
  objs = [u.ims, u.iris_ims, u.nlfff_ims]
  for i=0, n_elements(objs)-1 do begin
    objs[i].xr = u.xr
    objs[i].yr = u.yr
  endfor
  u.curx = x
  u.cury = y
  oWin.uvalue = u
  oWin.Refresh
  RETURN, 0
END

function mouseup, oWin, x, y, iButton, KeyMods, nClicks
  u = oWin.uvalue
  if ibutton eq 2 then u.buttondown = 0
  oWin.uvalue = u
  return, 0
end
      
;w03.uvalue = {ims:[im01, im02, im03], iris_ims:[im011, im021, im031], $
;              nlfff_ims:[im012, im022, im032], eis_xp:eis_xp, eis_yp0:eis_yp0, $
;              int_map:int_map, vel_map:vel_map, v_nt_map:v_nt_map, $
;              cos_i_iris:cos_i_iris, iris_xpix:iris_xpix, iris_ypix:iris_ypix, $
;              amp:amp, vel:vel, nth:nth, xpos:iris_xp, ypos:iris_yp, i:i, $
;              ion_list:ion_list, loc:loc, offset_y:offset_y, xr:xr, yr:yr}
  
function blink, window, $
  IsASCII, Character, KeyValue, X, Y, Press, Release, KeyMods
  
  window.refresh, /disable
  u = window.uvalue
  IF release THEN RETURN, 1
  ;  help, window, isascii, character, keyvalue, x, y, press, release
  character = string(character)
  if isASCII eq 1 then begin
    if total(character eq ['n', 'b', 'p']) then begin
      del = (character eq 'n') ? 1 : -1    
      i = (u.i + del) mod n_elements(u.int_map)
      u.i = i
      int_min = u.ims[0].min
      vel_min = u.ims[1].min
      vel_max = u.ims[1].max
      nth_min = u.ims[2].min
      nth_max = u.ims[2].max
      eis_yp = u.eis_yp0 - u.offset_y[i]
      u.ims[0].setdata, u.int_map[i].data, u.eis_xp, eis_yp
      u.ims[1].setdata, u.vel_map[i].data, u.eis_xp, eis_yp
      u.ims[2].setdata, u.v_nt_map[i].data, u.eis_xp, eis_yp
      hist = histogram(u.int_map[i].data, loc=xbin)
      pdf = total(hist, /cum)/total(finite(u.int_map[i].data))
      int_max = xbin[min(where(pdf gt 0.98))]
      u.ims[0].min = 0
      u.ims[0].max = int_max
      u.ims[1].min = vel_min
      u.ims[1].max = vel_max
      u.ims[2].min = nth_min
      u.ims[2].max = nth_max
      u.ims[0].title = 'Int: ' + u.ion_list[u.i].line_id + $
        ' (log T='+string(u.ion_list[u.i].logt, f='(f0.1)')+' K)'
      
    endif else if string(character) eq 'q' then begin
      u.ims[0].close
      return, 0
    endif else begin
      loc = (u.loc + 1) mod 3
      u.loc = loc
      if loc eq 0 then begin ; show nlfff
        for k=0, 2 do begin 
          u.ims[k].hide = 1
          u.iris_ims[k].hide = 1
          u.nlfff_ims[k].hide = 0
          u.nlfff_ims[k].order, /bring_to_front
        endfor
;      u.ims[0].title = 'Cos i (B vs. LOS)' 
      u.ims[0].title = 'B$_{LOS}$'
      endif else if loc eq 1 then begin ; show iris 
        for k=0, 2 do begin
          u.ims[k].hide = 1
          u.iris_ims[k].hide = 0
          u.nlfff_ims[k].hide = 1
          u.iris_ims[k].order, /bring_to_front
        endfor
      u.ims[0].title = 'Int: Si IV 1403 (log T=4.8 K)'        
      endif else if loc eq 2 then begin ; show eis
        for k=0, 2 do begin
          u.ims[k].hide = 0
          u.iris_ims[k].hide = 1
          u.nlfff_ims[k].hide = 1
          u.ims[k].order, /bring_to_front
        endfor
        u.ims[0].title = 'Int: ' + u.ion_list[u.i].line_id + $
          ' (log T='+string(u.ion_list[u.i].logt, f='(f0.1)')+' K)'
      endif
    endelse
  endif  
  window.uvalue = u
  window.refresh
  return, 0
end  
;; analysis
;

path = '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/20200405_120335'
restore, path+'/instr_info.sav'
if 0 then begin
  iris_path = (((iris_info.umodes).split('/'))[0:7]).join('/')
  iris_file = file_search(iris_path+'/*raster*.fits')
  if n_elements(iris_file) gt 1 then iris_file = iris_file[0]
  get_iris_si_iv_1403_fit, iris_file, amp, vel, nth, iris_xp, iris_yp
  save, iris_file, amp, vel, nth, iris_xp, iris_yp, $
        filename=path+'/iris_fit_res.sav'
endif else begin
  restore, path+'/iris_fit_res.sav'
endelse

restore, path+'/nlfff_results.sav'
read_sdo, blos_filename, dum, data_blos, /sil, /noshell, /use_shared_lib
get_xp_yp, index_blos[0], hmi_xp, hmi_yp, hmi_xxp, hmi_yyp, data=data_blos


iris_xpix = interpol(findgen(n_elements(res_xp)), res_xp, iris_xp)
iris_ypix = interpol(findgen(n_elements(res_yp)), res_yp, iris_yp)
;stop
cos_i_iris = interpolate(reform(abs(cos_i0[*, *, 6])), iris_xpix, iris_ypix, $
                         missing=!values.f_nan)

restore, path+'/EIS_fit_results.sav', /verbose
; int_map_all, vel_map_all, v_nt_map_all, ion_data_list
dum = where(ion_data_list.mass ne 0, /null)
order = sort(ion_data_list[dum].logt)
ion_list = (ion_data_list[dum])[order]
int_map = (int_map_all[dum])[order]
vel_map = (vel_map_all[dum])[order]
v_nt_map = (v_nt_map_all[dum])[order]
get_xp_yp, int_map[0], eis_xp, eis_yp0

i = 0
loc = 0
w03 = window(dim=[1e3, 5e2])
hist = histogram(int_map[i].data, loc=xbin)
pdf = total(hist, /cum)/total(finite(int_map[i].data))
int_max = xbin[min(where(pdf gt 0.99))]
xmar = 0.08
xs = (1.-6*xmar)/3.
xi = xmar + findgen(3)*(xs+2*xmar)
yi = 0.15
yf = 0.9

line_wave = !null
for j=0, n_elements(ion_list)-1 do begin
  wave_ = float((strsplit(ion_list[j].line_id, ' ', /extract))[-1])
  line_wave = [line_wave, wave_]
endfor
offset_y = eis_ccd_offset(line_wave)
eis_yp = eis_yp0 - offset_y[i]
xr = minmax(iris_xp)
yr = minmax(iris_yp)

im01 = IMAGE_(int_map[i].data, eis_xp, eis_yp, /current, aspect_r=1, $
              pos=[xi[0], yi, xi[0]+xs, yf], $
              TITLE='Int: ' + ion_list[i].line_id + $
              ' (log T='+string(ion_list[i].logt, f='(f0.1)')+' K)', $
              min=0, max=int_max, cb=cb01, xtickdir=1, ytickdir=1, $
              xr=minmax(iris_xp), yr=minmax(iris_yp), xtickinterval=50, $
              xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)')
im02 = IMAGE_(vel_map[i].data, eis_xp, eis_yp, /current, aspect_r=1, $
              pos=[xi[1], yi, xi[1]+xs, yf], $
              TITLE='Vel', $
              rgb_table=33, min=-30, max=30, cb=cb02, xtickdir=1, ytickdir=1, $
              xr=minmax(iris_xp), yr=minmax(iris_yp), xtickinterval=50, $
              xtitle='Solar X (arcsec)')
cb02.title = 'km/s'
im03 = IMAGE_(v_nt_map[i].data, eis_xp, eis_yp, /current, aspect_r=1, $
              pos=[xi[2], yi, xi[2]+xs, yf], $
              TITLE='V_nt', $
              rgb_table=4, min=0, max=60, cb=cb03, xtickdir=1, ytickdir=1, $
              xr=minmax(iris_xp), yr=minmax(iris_yp), xtickinterval=50, $
              xtitle='Solar X (arcsec)')
cb03.title = 'km/s'

im011 = image_(amp, iris_xp, iris_yp, over=im01, min=0, max=5e2)
im021 = image_(vel, iris_xp, iris_yp, over=im02, min=-3e1, max=3e1, rgb_table=33)
im031 = image_(nth, iris_xp, iris_yp, over=im03, min=0, max=3e1, rgb_table=4)

im013 = image_(cos_i0[*, *, 6], res_xp, res_yp, over=im01, min=0, max=1, rgb_table=22)
im023 = image_(cos_i0[*, *, 6], res_xp, res_yp, over=im02, min=0, max=1, rgb_table=22)
im033 = image_(cos_i0[*, *, 6], res_xp, res_yp, over=im03, min=0, max=1, rgb_table=22)
im013.title = 'Cos i (B vs. LOS)'

im012 = image_(data_blos, hmi_xp, hmi_yp, over=im01, min=-500, max=500, rgb_table=0)
im022 = image_(data_blos, hmi_xp, hmi_yp, over=im02, min=-500, max=500, rgb_table=0)
im032 = image_(data_blos, hmi_xp, hmi_yp, over=im03, min=-500, max=500, rgb_table=0)
im012.title = 'B$_{LOS}$'


w03.uvalue = {ims:[im01, im02, im03], iris_ims:[im011, im021, im031], $
              nlfff_ims:[im012, im022, im032], eis_xp:eis_xp, eis_yp0:eis_yp0, $
              int_map:int_map, vel_map:vel_map, v_nt_map:v_nt_map, $
              cos_i_iris:cos_i_iris, iris_xpix:iris_xpix, iris_ypix:iris_ypix, $
              amp:amp, vel:vel, nth:nth, xpos:iris_xp, ypos:iris_yp, i:i, $
              ion_list:ion_list, loc:loc, offset_y:offset_y, xr:xr, yr:yr, $
              buttondown:0L, curx:0l, cury:0l}
w03.keyboard_handler='blink'
w03.MOUSE_WHEEL_HANDLER = 'sync_fov_wheel'
w03.mouse_down_handler = 'mousedown'
w03.mouse_motion_handler = 'sync_fov_move'
w03.mouse_up_handler = 'mouseup'
                 
end
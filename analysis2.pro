
;w03.uvalue = {ims:[im01, im02, im03], iris_ims:[im011, im021, im031], $
;              nlfff_ims:[im012, im022, im032], eis_xp:eis_xp, eis_yp:eis_yp, $
;              int_map:int_map, vel_map:vel_map, v_nt_map:v_nt_map, $
;              cos_i_iris:cos_i_iris, iris_xpix:iris_xpix, iris_ypix:iris_ypix, $
;              amp:amp, vel:vel, nth:nth, xpos:xpos, ypos:ypos, i:i, $
;              ion_list:ion_list, loc:loc}
  
function blink, window, $
  IsASCII, Character, KeyValue, X, Y, Press, Release, KeyMods
  
  window.refresh, /disable
  u = window.uvalue
  IF release THEN RETURN, 1
  ;  help, window, isascii, character, keyvalue, x, y, press, release

  if isASCII eq 1 then begin
    if string(character) eq 'n' then begin 
      i = (u.i + 1) mod n_elements(u.int_map)
      u.ims[0].setdata, int_map[i].data, u.eis_xp, u.eis_yp
      u.ims[1].setdata, vel_map[i].data, u.eis_xp, u.eis_yp
      u.ims[2].setdata, v_nt_map[i].data, u.eis_xp, u.eis_yp
      u.ims[0].title = 'Int: ' + u.ion_list[i]
    endif else begin
      loc = (u.loc + 1) mod 3
      if loc eq 0 then begin ; show nlfff
        for k=0, 3 do begin 
          u.ims[k].hide = 1
          u.iris_ims[k].hide = 1
          u.nlfff_ims[k].hide = 0
        endfor
      endif else if loc eq 1 then begin ; show iris 
        for k=0, 3 do begin
          u.ims[k].hide = 1
          u.iris_ims[k].hide = 0
          u.nlfff_ims[k].hide = 1
        endfor
      endif else if loc eq 2 then begin ; show eis
        for k=0, 3 do begin
          u.ims[k].hide = 0
          u.iris_ims[k].hide = 1
          u.nlfff_ims[k].hide = 1
        endfor
      endif
    endelse
  u.i = i
  u.loc = loc
  window.uvalue = u
  window.refresh
  return, 0
end  
;; analysis
;

path = '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/20190412_145240'
restore, path+'/instr_info.sav'
if 0 then begin
  iris_path = (((iris_info.umodes).split('/'))[0:7]).join('/')
  iris_file = file_search(iris_path+'/*raster*.fits')
  get_iris_si_iv_1403_fit, iris_file, amp, vel, nth, xpos, ypos
  save, iris_file, amp, vel, nth, xpos, ypos, $
        filename=path+'/iris_fit_res.sav'
endif else begin
  restore, path+'/iris_fit_res.sav'
endelse

restore, path+'/nlfff_results.sav'
iris_xpix = interpol(findgen(n_elements(res_xp)), res_xp, xpos)
iris_ypix = interpol(findgen(n_elements(res_yp)), res_yp, ypos)
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
get_xp_yp, int_map[0], eis_xp, eis_yp

i = 0
loc = 0
w03 = window(dim=[1e3, 5e2])

im01 = IMAGE_(int_map[i].data, eis_xp, eis_yp, LAYOUT=[3,1,1], $
  TITLE='Int: '+ion_list[i].line_id, /current)
im02 = IMAGE_(vel_map[i].data, eis_xp, eis_yp, LAYOUT=[3,1,2], TITLE='Vel', $
  /CURRENT, rgb_table=33, min=-30, max=30)
im03 = IMAGE_(v_nt_map[i].data, eis_xp, eis_yp, LAYOUT=[3,1,3], TITLE='V_nt', $
  /CURRENT, rgb_table=4, min=0, max=30)

im011 = image_(amp, xpos, ypos, over=im01, min=0, max=5e2)
im021 = image_(vel, xpos, ypos, over=im02, min=-3e2, max=3e2, rgb_table=33)
im031 = image_(nth, xpos, ypos, over=im03, min=0, max=5e2, rgb_table=4)

im012 = image_(cos_i_iris, xpos, ypos, over=im01, min=0, max=1, rgb_table=22)
im022 = image_(cos_i_iris, xpos, ypos, over=im02, min=0, max=1, rgb_table=22)
im032 = image_(cos_i_iris, xpos, ypos, over=im03, min=0, max=1, rgb_table=22)

w03.uvalue = {ims:[im01, im02, im03], iris_ims:[im011, im021, im031], $
              nlfff_ims:[im012, im022, im032], eis_xp:eis_xp, eis_yp:eis_yp, $
              int_map:int_map, vel_map:vel_map, v_nt_map:v_nt_map, $
              cos_i_iris:cos_i_iris, iris_xpix:iris_xpix, iris_ypix:iris_ypix, $
              amp:amp, vel:vel, nth:nth, xpos:xpos, ypos:ypos, i:i, $
              ion_list:ion_list, loc:loc}
w03.keyboard_handler='blink'

                 
end
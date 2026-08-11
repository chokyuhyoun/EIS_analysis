pro nlfff_cal, cen, xr, yr, ref_time_str0, save_current_path = save_current_path

cd, current=path
if ~keyword_set(save_current_path) then begin
  save_path = path + path_sep() + ref_time_str0
  file_mkdir, save_path
  cd, save_path
endif else begin
  save_path = path      
endelse


; Select the FOV

xbi = cen[0] - 0.5*xr
xbf = cen[0] + 0.5*xr
ybi = cen[1] - 0.5*yr
ybf = cen[1] + 0.5*yr
xbcen = cen[0]
ybcen = cen[1]


npix_grid = 50 ; just for drawing 
equi_dist = 725.*0.5 ; = HMI a pixel size at the disk center. 725 km/arcsec * 0.5 arcsec/pix
;ref_time_str = '2020-02-02 02:02:02' ; reference time
ref_time_str = anytim2utc(ref_time_str0, /ccsds)
ref_time = anytim(ref_time_str) 

; Download/read data. 
; If you use the LMSAL ethernet, just read them directly from JSOC server via /files_only keyword.
; If not, remove /files_only keyword and download them. 
files_only = 1
 
ssw_jsoc_time2data, anytim(ref_time-360, /ccsds), anytim(ref_time+360, /ccsds), $
                    index_i, ds='hmi.B_720s', segment='inclination', $
                    files_only=files_only, inclination_file, /silent
ssw_jsoc_time2data, anytim(ref_time-360, /ccsds), anytim(ref_time+360, /ccsds), $
                    index_a, ds='hmi.B_720s', segment='azimuth', $
                    files_only=files_only, azimuth_file, /silent
ssw_jsoc_time2data, anytim(ref_time-360, /ccsds), anytim(ref_time+360, /ccsds), $
                    index_b, ds='hmi.B_720s', segment='field', $
                    files_only=files_only, field_file, /silent
ssw_jsoc_time2data, anytim(ref_time-360, /ccsds), anytim(ref_time+360, /ccsds), $
                    index_d, ds='hmi.B_720s', segment='disambig', $
                    files_only=files_only, disambig_file, /silent
ssw_jsoc_time2data, anytim(ref_time-6, /ccsds), anytim(ref_time+6, /ccsds), $
                    aia_index, ds='aia.lev1_euv_12s', wave=171, $
                    files_only=files_only, aia171_file, /silent

;azimuth_file = file_search(hmi_path, '*azimuth.fits')
;inclination_file = file_search(hmi_path, '*inclination.fits')
;field_file = file_search(hmi_path, '*field.fits')
;disambig_file = file_search(hmi_path, '*disambig.fits')
read_sdo, azimuth_file[0], dum, azi_data, /sil, /noshell, /use_shared_lib
read_sdo, inclination_file[0], dum, inc_data, /sil, /noshell, /use_shared_lib
read_sdo, field_file[0], dum, field_data, /sil, /noshell, /use_shared_lib
read_sdo, disambig_file[0], dum, disamb_data, /sil, /noshell, /use_shared_lib
read_sdo, aia171_file[0], dum, aia171_data, /sil, /noshell, /use_shared_lib

hmi_time = anytim2tai(index_i.date_obs)

print, systime()
init_time = systime(/sec)
inc_data = rotate(temporary(inc_data), 2)
azi_data = rotate(temporary(azi_data), 2)
field_data = rotate(temporary(field_data), 2)
disamb_data = rotate(temporary(disamb_data), 2)
hmi_disambig, azi_data, disamb_data, 2

aia_xp = (findgen(aia_index.naxis1) - aia_index.crpix1 + 1)*aia_index.cdelt1 + aia_index.crval1
aia_yp = (findgen(aia_index.naxis2) - aia_index.crpix2 + 1)*aia_index.cdelt2 + aia_index.crval2

index = index_b
index_i0 = index_i
index.crpix1 = index.naxis1 - index.crpix1 + 1
index.crpix2 = index.naxis2 - index.crpix2 + 1
r_sun = index.rsun_obs

  ; x : toward the observer / y : solar west /  z : solar north  
bx0 = field_data*cos(inc_data*!dtor)
by0 = field_data*sin(inc_data*!dtor)*sin(azi_data*!dtor)
bz0 = -field_data*sin(inc_data*!dtor)*cos(azi_data*!dtor)
  
hmi_xp = (findgen(index.naxis1) - index.crpix1 + 1)*index.cdelt1 + index.crval1
hmi_yp = (findgen(index.naxis2) - index.crpix2 + 1)*index.cdelt2 + index.crval2
hmi_xxp = rebin(hmi_xp, index.naxis1, index.naxis2)
hmi_yyp = rebin(transpose(hmi_yp), index.naxis1, index.naxis2)

from_rect = transpose([[(sqrt(r_sun^2.-hmi_xxp^2.-hmi_yyp^2.))[*]], [hmi_xxp[*]], [hmi_yyp[*]]])
hmi_sphere = cv_coord(from_rect=from_rect, /to_sphere)
  
hmi_lon = reform(hmi_sphere[0, *], index.naxis1, index.naxis2)
hmi_lat = reform(hmi_sphere[1, *], index.naxis1, index.naxis2)
 
br0 = bx0*cos(hmi_lat)*cos(hmi_lon) + by0*cos(hmi_lat)*sin(hmi_lon) + bz0*sin(hmi_lat)
blon0 = -bx0*sin(hmi_lon) + by0*cos(hmi_lon)
blat0 = -bx0*sin(hmi_lat)*cos(hmi_lon) - by0*sin(hmi_lat)*sin(hmi_lon) + bz0*cos(hmi_lat)

; check boundary position in lon, lat
dum_y = [xbi, xbf, xbcen, xbcen]
dum_z = [ybcen, ybcen, ybf, ybi]
dum_x = sqrt(r_sun^2. - dum_y^2. - dum_z^2.)
from_rect1 = transpose([[dum_x], [dum_y], [dum_z]]) 
dum_res = cv_coord(from_rect=from_rect1, /to_sphere) ; in radian
hel_left = dum_res[*, 0]
hel_right = dum_res[*, 1]
hel_top = dum_res[*, 2]
hel_bottom = dum_res[*, 3]
hel_cenx = 0.5*(hel_left[0] + hel_right[0])
hel_ceny = 0.5*(hel_top[1] + hel_bottom[1])
  
; equi degree grid
pix2latrad = equi_dist/(r_sun*725.)
nlat_pix = floor(round((hel_top[1] - hel_bottom[1])/pix2latrad)/npix_grid)*npix_grid + 1
pix2lonrad = equi_dist/(r_sun*cos(hel_ceny)*725.)
nlon_pix = floor(round((hel_right[0] - hel_left[0])/pix2lonrad)/npix_grid)*npix_grid + 1

latp = (findgen(nlat_pix) - 0.5*(nlat_pix - 1))*pix2latrad + hel_ceny
latpp = rebin(transpose(latp), nlon_pix, nlat_pix)

lonpp = fltarr(nlon_pix, nlat_pix)
for dum = 0, nlat_pix-1 do $
  lonpp[*, dum] = (findgen(nlon_pix) - 0.5*(nlon_pix-1)) $
                  *equi_dist/(r_sun*cos(latp[dum])*725.) + hel_cenx 
equi_dist_co0 = cv_coord(from_sphere=transpose([[lonpp[*]], [latpp[*]], [fltarr(n_elements(lonpp))+r_sun]]), $
                         /to_rect)
equi_dist_co = reform(equi_dist_co0[1:*, *], 2, nlon_pix, nlat_pix)
  
interp_xp = interpol(findgen(index.naxis1), hmi_xp, reform(equi_dist_co[0, *, *]))
interp_yp = interpol(findgen(index.naxis2), hmi_yp, reform(equi_dist_co[1, *, *])) 
br0_equi_dist = interpolate(br0, interp_xp, interp_yp)
blon0_equi_dist = interpolate(blon0, interp_xp, interp_yp)
blat0_equi_dist = interpolate(blat0, interp_xp, interp_yp)

save, br0_equi_dist, blon0_equi_dist, blat0_equi_dist, index_i0, filename=save_path+'/nlfff_input.sav'

res_file = file_search(save_path, 'field.dat', count=n)
if n eq 0 then begin
  optimization_fff, blon0_equi_dist, blat0_equi_dist, br0_equi_dist
endif
restore, save_path+'/nlfff_input.sav'
Read_bfield_fff, 'field.dat', xsize, ysize, zsize, bx, by, bz, x, y, z  ; lon, lat, r
norm_size = max([n_elements(x), n_elements(y), n_elements(z)])
z = findgen(n_elements(z))

; drawing the field lines
dfp = 30
fpx = [0.5*dfp : nlon_pix-0.5*dfp : dfp]
fpy = [0.5*dfp : nlat_pix-0.5*dfp : dfp]
fpxx = rebin(fpx, n_elements(fpx), n_elements(fpy))
fpyy = rebin(transpose(fpy), n_elements(fpx), n_elements(fpy))

real = where(abs(br0_equi_dist[fpxx[*], fpyy[*]]) gt 100)
fpxx = fpxx[real]
fpyy = fpyy[real]
lines = list(!null)
for dum = 0, n_elements(fpxx)-1 do begin &$    
  b_line, bx, by, bz, [fpxx[dum], fpyy[dum], 0], r &$ ; in normalized coordinate
  lines.add, r &$
endfor 
lines.remove, 0

lines_lonlat = list(!null)
for dum = 0, n_elements(lines)-1 do begin &$
  dum_lon = interpolate(lonpp, (lines[dum])[*, 0], (lines[dum])[*, 1]) &$
  dum_lat = interpol(latp, findgen(nlat_pix), (lines[dum])[*, 1]) &$
  dum_r = interpol(z, findgen(n_elements(z)), (lines[dum])[*, 2]) &$
  res = cv_coord(from_sphere=transpose([[dum_lon], [dum_lat], [dum_r+r_sun]]), /to_rect) &$
  lines_lonlat.add, res &$
endfor
lines_lonlat.remove, 0

n_h = 20
; grid within xbi-xbf, ybi-ybf with equi dist(degree)
res_xp = hmi_xp[where((hmi_xp gt min(equi_dist_co[0, *, *])) and $
                      (hmi_xp lt max(equi_dist_co[0, *, *])), /null)] 
res_yp = hmi_yp[where((hmi_yp gt min(equi_dist_co[1, *, *])) and $
                      (hmi_yp lt max(equi_dist_co[1, *, *])), /null)]
res_zp = findgen(n_h)*index.cdelt1   ;should be radial? or los?
res_xxp = rebin(res_xp, n_elements(res_xp), n_elements(res_yp), n_h)
res_yyp = rebin(transpose(res_yp), n_elements(res_xp), n_elements(res_yp), n_h)
res_zzp0 = rebin(reform(res_zp, 1, 1, n_h), n_elements(res_xp), n_elements(res_yp), n_h)
res_zzp = sqrt((r_sun+res_zzp0)^2. - res_xxp^2. -res_yyp^2.)

from_rect2 = transpose([[res_zzp[*]], [res_xxp[*]], [res_yyp[*]]]) 
dum_res2 = cv_coord(from_rect=from_rect2, /to_sphere) ; [lon, lat, r]
res2_lon = reform(dum_res2[0, *], n_elements(res_xp), n_elements(res_yp), n_h)
res2_lat = reform(dum_res2[1, *], n_elements(res_xp), n_elements(res_yp), n_h)
res2_latpp = (res2_lat - hel_ceny)/pix2latrad + 0.5*(nlat_pix-1)
res2_lonpp = res2_lon*0.
for dumz = 0, n_h-1 do for dum = 0, n_elements(res_yp)-1 do $
    res2_lonpp[*, dum, dumz] = (res2_lon[*, dum, dumz] - hel_cenx)/(equi_dist/725./(cos(res2_lat[0, dum, dumz])*(r_sun+res_zp[dumz]))) + 0.5*(nlon_pix-1)

cos_i0 = fltarr(n_elements(res_xp), n_elements(res_yp), n_h)
b_obs0 = cos_i0
for ii=0, n_h-1 do begin 
  bx_obs = bz[*, *, ii]*cos(latpp)*cos(lonpp) $ ; los
          - bx[*, *, ii]*sin(lonpp) $
          - by[*, *, ii]*sin(latpp)*cos(lonpp)
  by_obs = bz[*, *, ii]*cos(latpp)*sin(lonpp) $ ; west
          + bx[*, *, ii]*cos(lonpp) $
          - by[*, *, ii]*sin(latpp)*sin(lonpp)
  bz_obs = bz[*, *, ii]*sin(latpp) $            ; north
          + by[*, *, ii]*cos(latpp)
  b_obs = sqrt(bx_obs^2. + by_obs^2. + bz_obs^2.)
  cos_i = bx_obs/b_obs
  br_sgn = sgn(interpolate(bz[*, *, ii], res2_lonpp[*, *, ii], res2_latpp[*, *, ii], missing=!values.f_nan))
  cos_i0[*, *, ii] = interpolate(cos_i, res2_lonpp[*, *, ii], res2_latpp[*, *, ii], missing=!values.f_nan)*br_sgn
  ; now cos_i0 > 0 --> i0 < 90 degree, cos_i0 < 0 --> i0 > 90 degree
  b_obs0[*, *, ii] = interpolate(b_obs, res2_lonpp[*, *, ii], res2_latpp[*, *, ii], missing=!values.f_nan)

endfor
save, cos_i0, res_xp, res_yp, index_i0, filename=save_path+'/nlfff_results.sav'

w01 = window(dim=[1d3, 1d3], /buffer)
im01 = image(bx0, hmi_xp, hmi_yp, /current, axis=2, pos=[0.1, 0.6, 0.4, 0.9], $
              min=-1d3, max=1d3, title='(a) B LOS and equidistant grid', $
              xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
              xr=xbcen+200.*[-1, 1], yr=ybcen+200.*[-1, 1])
t01 = text(im01.pos[0]+0.02, im01.pos[3]-0.02, strmid(index_i0.date_obs.replace('T', ' '), 0, 19)+' UT', $
           /normal, color='black', font_size=13, align=0, vertical_align=1)
for dum0 = 0, nlon_pix-1, npix_grid do $
  p011 = plot(equi_dist_co[0, dum0, *], equi_dist_co[1, dum0, *], 'y', over=im01, transp=50)
for dum0 = 0, nlat_pix-1, npix_grid do $
  p021 = plot(equi_dist_co[0, *, dum0], equi_dist_co[1, *, dum0], 'y', over=im01, transp=50)

im02 = plot3d(/test, /nodata, /current, axis_style=2, $
              aspect_ratio=1, aspect_z=1, $
              ;title='B$_r$ and NLFFF lines', $
              xtitle='Solar East (pix)', ytitle='Solar North (pix)', $
              xr=[0, nlon_pix], yr=[0, nlat_pix], zr=[0, n_elements(z)])
for ii=2, 11 do im02.axes[ii].hide=1
im02.axes[8].hide = 0
im02.axes[8].title = 'Normal to Solar Surface (pix)'
im02.axes[8].tickfont_size = 10
im02.axes[8].showtext = 1
im02.axes[8].minor = 1
im02.axes[1].minor = 1
im02.axes[1].tickinterval = 200
im02.axes[8].tickinterval = 200
  
im02.scale, 0.3, 0.3, 0.3
im02.pos = im02.pos + [0.8, 0.85, 0.8, 0.85]              
t02 = text(0.72, 0.9, '(b) B$_r$ and NLFFF lines', align=0.5, /normal, font_size=13)
im021 = image(br0_equi_dist, zvalue=0, min=-1d3, max=1d3, over=im02)                
for dum = 0, n_elements(fpxx)-1 do $
    im022 = plot3d([(lines[dum])[*, 0]], [(lines[dum])[*, 1]], [(lines[dum])[*, 2]], over=im02, $
                color='g', thick=1, transp=0)
  
im03 = image(aia171_data, aia_xp, aia_yp, /current, axis=2, pos=[0.1, 0.1, 0.4, 0.4], $
              min=0, max=1.5d3, rgb_table=aia_lct_fun(171), $
              xr=xbcen+200.*[-1, 1], yr=ybcen+200.*[-1, 1], $
              xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
              title='(c) AIA 171 $\AA$ and NLFFF lines')
t03 = text(im03.pos[0]+0.02, im03.pos[3]-0.02, strmid(aia_index.date_obs.replace('T', ' '), 0, 19)+' UT', $
                /normal, color='white', font_size=13, align=0, vertical_align=1)
for dum = 0, n_elements(lines_lonlat)-1 do $
  p1240 = plot((lines_lonlat[dum])[1, *], (lines_lonlat[dum])[2, *], over=im03, 'g', transp=0)

height_pix = 6
im125 = image(cos_i0[*, *, height_pix], res_xp, res_yp, /current, axis=2, $
              pos=[0.6, 0.1, 0.9, 0.4], $
               xr=xbcen+200.*[-1, 1], yr=ybcen+200.*[-1, 1], $
               title='(d) Cos !8i!x = '+string(height_pix*0.5*725., f='(i0)')+' km', $
               xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
               rgb_table=20, min=0, max=1)
cb125 = colorbar(target=im125, orientation=1, border=1, textpos=1, $
                  pos = im125.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0], title='Cos !8i!x')       
;w01.save, save_path+'/NLFFF_result.png', resol=200, /bitmap     
w01.save, save_path+'/NLFFF_result.pdf', page_size=w01.dimen/1d2, width=w01.dimen[1]/1d2, /bitmap
w01.close       
save, cos_i0, height_pix, res_xp, res_yp, filename=save_path+'/cos_i0.sav'
cd, path
print, string((systime(/sec) - init_time)/60., f='(f0.2)')+' min'
end

;time_str = ['20140704_114000', $
;            '20140705_230030', $
;            '20150224_190314', $
;            '20190410_121535', $
;            '20190412_145240', $
;            '20190413_013438', $
;            '20190413_113910', $
;            '20190417_011605', $
;            '20190417_154227', $
;            '20200402_224709', $
;            '20200403_011735', $
;            '20200404_021153', $
;            '20200405_120335', $
;            '20200405_184659', $
;            '20200406_164019', $
;            '20200406_233539', $
;            '20200407_025456', $
;            '20190415_124334', $
;            '20190415_215330']
;center = [[-107.5, -243.5], $
;          [139.5, -244.5], $
;          [37.0, -40.0], $
;          [-685.5, 174.0], $
;          [-310.0, 197.5], $
;          [-227.5, 200.5], $
;          [-140.5, 202.0], $
;          [595.0, 168.0], $
;          [684.0, 158.5], $
;          [-401.0, 540.0], $
;          [-383.0, 542.0], $
;          [-199.0, 547.5], $
;          [68.0, 548.5], $
;          [124.5, 548.5], $
;          [289.5, 546.0], $
;          [344.5, 540.0], $
;          [363.0, 537.5], $
;          [324.5, 184.5], $
;          [392.0, 183.5]] 
;xr = [269, 399, 440, 331, 440, 453, 459, 388, 332, 314, $
;      312, 318, 296, 283, 249, 235, 232, 493, 468]*1.              
;yr = [209, 215, 278, 300, 281, 275, 270, 236, 231, 174, $
;      170, 155, 157, 157, 158, 144, 145, 251, 247]*1.  
;
;;for i=0, n_elements(xr)-1 do begin
;for i=2, 2 do begin
;  print, time_str[i], center[*, i], xr[i], yr[i]
;  nlfff_cal, center[*, i], xr[i], yr[i], time_str[i]
;endfor
nlfff_cal, [-275., -75.], 175., 225., '2013-10-07 02:10:00', /save_current_path
end
      
      
      

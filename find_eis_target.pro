eis_t0 = '2013-07-16T00:00:00'
eis_t1 = '2026-07-28T00:00:00

if 0 then begin 
  str = eis_obs_structure(eis_t0, eis_t1, /quick, /quiet)
  save, str, filename='/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/eis_str.sav'
endif 
restore, '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/eis_str.sav'  
c_slit = str.slit_width le 2
c_exptime = float(str.exp_times) gt 10.
c_sci_obj = strmatch(str.sci_obj, '*AR*')
c_fovx = str.fovx gt 50.
c_r_sun = sqrt(str.xcen^2. + str.ycen^2.) lt 600.
c_nwin = ((str.wavelnth.strlen()+1.)/6.) gt 10

c_all = where(c_slit and c_exptime and c_sci_obj and c_fovx and c_r_sun and c_nwin, /null)
eis_str = str[c_all]


w01 = window(dim=[8d2, 8d2])
aia_lct, rr, gg, bb, wavelnth=171
ct171 = [[rr], [gg], [bb]]

for i = 6012, n_elements(eis_str)-1 do begin
  print, i, n_elements(eis_str), i*1d2/n_elements(eis_str), $
        f='(i6, "/", i6, "   ", f5.1, " %")'
  eis_obs_t0 = eis_str[i].date_obs
  eis_obs_t1 = eis_str[i].date_end
  eis_obs_tm = anytim(0.5*(anytim(eis_obs_t0) + anytim(eis_obs_t1)), /ccsds)
  eis_x0 = eis_str[i].xcen - eis_str[i].fovx*0.5
  eis_x1 = eis_str[i].xcen + eis_str[i].fovx*0.5
  eis_y0 = eis_str[i].ycen - eis_str[i].fovy*0.5
  eis_y1 = eis_str[i].ycen + eis_str[i].fovy*0.5
  
  iris_t0 = anytim(reltime(eis_obs_t0, hour=-6), /ccsds)
  iris_t1 = anytim(reltime(eis_obs_t1, hour=6), /ccsds)
  query = 'https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&'+ $
    'startTime='+iris_t0+'&stopTime='+iris_t1+'&hasData=true&hideMostLimbScans=true&limit=200'
  iris_info = ssw_hcr_query(query)
  
  d_all = !null
  for j = 0, n_elements(iris_info)-1 do begin
    if total(strmatch(tag_names(iris_info[j]), 'RASTER_FOVX')) eq 0 then continue
    d_fovx = iris_info[j].raster_fovx gt 30.
    d_t_overlap = (anytim(iris_info[j].starttime) lt anytim(eis_obs_t1)) and $
                  (anytim(iris_info[j].stoptime) gt anytim(eis_obs_t0))
  
    iris_x0 = iris_info[j].xcen - iris_info[j].raster_fovx*0.5
    iris_x1 = iris_info[j].xcen + iris_info[j].raster_fovx*0.5
    iris_y0 = iris_info[j].ycen - iris_info[j].raster_fovy*0.5
    iris_y1 = iris_info[j].ycen + iris_info[j].raster_fovy*0.5
    d_x_overlap = (iris_x0 lt eis_x1) and (iris_x1 gt eis_x0)
    d_y_overlap = (iris_y0 lt eis_y1) and (iris_y1 gt eis_y0)
;    print, d_fovx, d_t_overlap,  d_x_overlap,  d_y_overlap
    d_all = [d_all, (d_fovx and d_t_overlap and d_x_overlap and d_y_overlap)] 
  endfor
  if n_elements(d_all) eq 0 then d_all = 0
  if total(d_all) gt 0 then begin
    iris_str = iris_info[where(d_all, /null)]
    ssw_jsoc_time2data, eis_obs_tm, reltime(eis_obs_tm, min=1), $
                        index_aia0, ds='aia.lev1_euv_12s', wave=171, $
                        files_only=1, filenames, /silent
    if n_elements(index_aia0) eq 0 then break
    index_aia0 = index_aia0[0]
    read_sdo, filenames[0], dum, data_aia0, /sil, /noshell, /use_shared_lib
    get_xp_yp, index_aia0, xp_aia0, yp_aia0
    eis_fov = eis_str[i].fovx > eis_str[i].fovy
    im01 = image_(data_aia0, xp_aia0, yp_aia0, $
                  pos = [0.1, 0.14, 0.9, 0.94], $
                  xr=eis_str[i].xcen + 0.5*(eis_fov+100)*[-1, 1], $
                  yr=eis_str[i].ycen + 0.5*(eis_fov+100)*[-1, 1], $
                  min=0, max=2e3, /current, rgb_table=ct171, $
                  title='AIA 171 $\AA$ '+index_aia0.date_obs)
    for k=0, n_elements(iris_str)-1 do begin
      p02 = plot(iris_str[k].xcen+0.5*iris_str[k].xfov*[-1, -1, 1, 1, -1], $
                 iris_str[k].ycen+0.5*iris_str[k].yfov*[-1, 1, 1, -1, -1], over=im01, $
                 'g2', trans=0)
      t02 = text(0.05, 0.06-k*0.03, 'IRIS: '+iris_str[k].starttime+' - '+iris_str[k].stoptime, $
                 color='green', font_size=13)
    endfor
  
    p03 = plot(eis_str[i].xcen + 0.5*eis_str[i].fovx*[-1, -1, 1, 1, -1], $
               eis_str[i].ycen + 0.5*eis_str[i].fovy*[-1, 1, 1, -1, -1], over=im01, $
               'c2', trans=0)
    t03 = text(0.05, 0.09, 'EIS: '+eis_str[i].date_obs+' - '+eis_str[i].date_end, $
               color='cyan', font_size=13)
    w01.save, '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff/all_targets/'+$
              string(i, f='(i04)')+'.png', resol=150
;      stop
  endif               
;  stop
  w01.erase
endfor               
end
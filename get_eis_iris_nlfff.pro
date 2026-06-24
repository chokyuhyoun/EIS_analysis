nlfff_path = '/System/Volumes/Data/links/mimas/sanhome2/khcho/Andy_nlfff'
f = file_search(nlfff_path, 'nlfff_results.sav')

;for i=0, n_elements(f)-1 do begin
  i = 14
  cd, file_dirname(f[i])
  print, f[i]
  restore, f[i]

  w01 = window(dim=[8e2, 8e2])
  t0 = index_i0.date_obs
  t1 = anytim(reltime(index_i0.date_obs, min=10), /ccsds)

  if 1 then begin
    ssw_jsoc_time2data, t0, t1, $
      index_blos, ds='hmi.M_45s', segment='magnetogram', $
      files_only=1, filenames, /silent
    index_blos = index_blos[0]
    read_sdo, filenames[0], dum, data_blos, /sil, /noshell, /use_shared_lib
  endif else begin
    filenames = file_search(file_dirname(f[i])+'/HMI*.fits')
    read_sdo, filenames[0], index_blos, data_blos, /sil 
  endelse
  
  get_xp_yp, index_blos[0], hmi_xp, hmi_yp, hmi_xxp, hmi_yyp, data=data_blos
  nlfff_xr = res_xp[-1] - res_xp[0]
  nlfff_yr = res_yp[-1] - res_yp[0]
  nlfff_fov = max([nlfff_xr, nlfff_yr])
  im01 = image_(data_blos, hmi_xp, hmi_yp, /current, $
                xr=mean(res_xp)+(nlfff_fov*0.5+50)*[-1, 1], $
                yr=mean(res_yp)+(nlfff_fov*0.5+50)*[-1, 1], $
                min=-500, max=500, xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)')
  p01 = plot((minmax(res_xp))[[0, 1, 1, 0, 0]], $  
             (minmax(res_yp))[[0, 0, 1, 1, 0]], over=im01, $
             'r2', trans=50)
  t01 = text(0.05, 0.1, 'NLFFF: '+index_i0.date_obs, color='red', font_size=13)
  query = 'https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&'+ $
    'startTime='+t0+'&stopTime='+t1+'&hasData=true&hideMostLimbScans=true&limit=200'
  iris_info = (ssw_hcr_query(query))[0]
  p02 = plot(iris_info.xcen+0.5*iris_info.xfov*[-1, -1, 1, 1, -1], $
             iris_info.ycen+0.5*iris_info.yfov*[-1, 1, 1, -1, -1], over=im01, $
             'g2', trans=50)
  t02 = text(0.05, 0.08, 'IRIS: '+iris_info.starttime, color='green', font_size=13)

  str = eis_obs_structure(reltime(t0,hours=-12), reltime(t0,hours=12), /quick, /quiet)
  dist_from_iris = sqrt((str.xcen - iris_info.xcen)^2. + (str.ycen - iris_info.ycen)^2.)
  eis_info = str[where(dist_from_iris eq min(dist_from_iris), /null)]
  eis_info = str[-5]
  p03 = plot(eis_info.xcen + 0.5*eis_info.fovx*[-1, -1, 1, 1, -1], $
             eis_info.ycen + 0.5*eis_info.fovy*[-1, 1, 1, -1, -1], over=im01, $
             'b2', trans=50)
  t03 = text(0.05, 0.06, 'EIS: '+eis_info.date_obs, color='blue', font_size=13)
  
;  md = vso_search(eis_info.date_obs, eis_info.date_end, instr='eis')
;  s = vso_get(md)
  w01.save, 'FOV_comp.png', resol=200
  w01.close
  blos_filename = filenames[0]
  save, iris_info, eis_info, index_blos, blos_filename, filename='instr_info.sav'
;stop
;endfor


end
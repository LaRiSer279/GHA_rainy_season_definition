load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "./functions.ncl"

begin
  id = "./" 
  in = addfile(id + "Filter_CHIRPS.nc","r")
	
;;;;; ==================== Read in variables ============================= ;;;;;;;;;
  rainfall = in->filtered_rainfall({-12.5:20.},{21.5:52},:,:) ; Rainfall data
  lon = in->lon({21.5:52})
  lat = in->lat({-12.5:20.})
  ref_mon = in->ref_month({-12.5:20.},{21.5:52})

;;;;; ==================== Pre-initialization ============================ ;;;;;;;;;;;;
rclim = new((/dimsizes(lat), dimsizes(lon), 3, 3, 2, 2/), "integer", -999) ;SAVE
rclim@units = "days since 1990-01-01"
rclim!0 = "lat"
rclim&lat = lat
rclim!1 = "lon"
rclim&lon = lon

rclim1 = rclim	;SAVE
rclim2 = rclim	;SAVE

copy_VarMeta(rclim, rclim1)
copy_VarMeta(rclim, rclim2)
 
thres = new((/dimsizes(lat), dimsizes(lon), 3/), "float", -999)	;SAVE
copy_VarCoords(rclim,thres)

;;;; ====================================================================== ;;;;;;;;;;;;; 
dimLon = dimsizes(lon)-1
dimLat = dimsizes(lat)-1

ind2 = ind_nearest_coord(21.125, lon, 0)
  ind1 = ind_nearest_coord(19.875, lat, 0) 
do ll = 0,dimLon
	do la = 0,dimLat
		print("Lon:" + lon(ll))
		print("Lat:" + lat(la))
		print(ref_mon(la,ll))
		
		;;;;; Skip ocean grid points ;;;;;;;;
		if(ismissing(ref_mon(la,ll)))then
			continue
		end if	
		
		; pentad rainfall climatology
		RR_new = dim_avg_n_Wrap(rainfall(la,ll,:,:), 0)	
		;print(RR_new)
	
		; define threshold for rainy season
		bb = get_threshold(RR_new, True)
		bb50 = bb[0]
		print(bb50)
		
		; define rainy season
		onset_clim50 = climatology_intervals(RR_new, bb50)
		print(onset_clim50)
		
		
		; confidence interval edges
		RR_perc = bootstrapping(rainfall(la,ll,:,:))
		bb5x = RR_perc(0,:)
		bb95x = RR_perc(1,:)

		bb5 = bb5x(bb[1])
		bb95 = bb95x(bb[1])
		print(bb5)
		print(bb95)
		onset_clim5 = climatology_intervals(RR_new, bb5)
		print(onset_clim5)
		onset_clim95 = climatology_intervals(RR_new, bb95)
		print(onset_clim95)
		
		thres(la,ll,:) = (/bb5, bb50, bb95/)


		ref_date = toint(cd_inv_calendar(1981,ref_mon(la,ll),1,0,0,0,"days since 1981-01-01",0))
		
		if(any(.not. ismissing(onset_clim5)))then
			if(dimsizes(dimsizes(onset_clim5)) .eq. 1)then
				rclim(la,ll,0,0,0,:) = onset_clim5 
				rclim(la,ll,0,0,1,:) = ref_date + 5*onset_clim5	
			else
				nseas = dimsizes(onset_clim5(:,0))-1
				onset = ref_date + 5*onset_clim5
			
				do seas = 0,nseas
					if(nseas .eq. 1)then
						rclim1(la,ll,0,seas,0,:) = (/onset_clim5(seas,:)/)
						rclim1(la,ll,0,seas,1,:) = (/onset(seas,:)/)
					else
						if(nseas .eq. 2)then
							rclim2(la,ll,0,seas,0,:) = (/onset_clim5(seas,:)/)
							rclim2(la,ll,0,seas,1,:) = (/onset(seas,:)/)
						end if
					end if		
				end do	
			
				delete(nseas)
				delete(onset)
			end if	
			delete(onset_clim5)
			delete(bb5)
		end if
		
		if(any(.not. ismissing(onset_clim50)))then	
			if(dimsizes(dimsizes(onset_clim50)) .eq. 1)then
				rclim(la,ll,1,0,0,:) = onset_clim50 
				rclim(la,ll,1,0,1,:) = ref_date + 5*onset_clim50 	
			else
				nseas = dimsizes(onset_clim50(:,0))-1
				onset = ref_date + 5*onset_clim50
			
				do seas = 0,nseas
					if(nseas .eq. 1)then
						rclim1(la,ll,1,seas,0,:) = (/onset_clim50(seas,:)/)
						rclim1(la,ll,1,seas,1,:) = (/onset(seas,:)/)
					else
						if(nseas .eq. 2)then
							rclim2(la,ll,1,seas,0,:) = (/onset_clim50(seas,:)/)
							rclim2(la,ll,1,seas,1,:) = (/onset(seas,:)/)
						end if
					end if		
				end do	
			
				delete(nseas)
				delete(onset)
			end if	
			delete(onset_clim50)
			delete(bb50)
		end if
			
		if(any(.not. ismissing(onset_clim95)))then
			if(dimsizes(dimsizes(onset_clim95)) .eq. 1)then
				rclim(la,ll,2,0,0,:) = onset_clim95 
				rclim(la,ll,2,0,1,:) = ref_date + 5*onset_clim95 	
			else
				nseas = dimsizes(onset_clim95(:,0))-1
				onset = ref_date + 5*onset_clim95
			
				do seas = 0,nseas
					if(nseas .eq. 1)then
						rclim1(la,ll,2,seas,0,:) = (/onset_clim95(seas,:)/)
						rclim1(la,ll,2,seas,1,:) = (/onset(seas,:)/)
					else
						if(nseas .eq. 2)then
							rclim2(la,ll,2,seas,0,:) = (/onset_clim95(seas,:)/)
							rclim2(la,ll,2,seas,1,:) = (/onset(seas,:)/)
						end if
					end if		
				end do	
			
				delete(nseas)
				delete(onset)
			end if	
			delete(onset_clim95)
			delete(bb95)
		end if
		
		delete(RR_new)
	end do
end do

tmp = rclim(:,:,:,:,1,:)
tmp_uni = where(tmp .gt. 365, tmp-365, tmp)
rclim(:,:,:,:,1,:) = (/tmp_uni/)
tmp = rclim1(:,:,:,:,1,:)
tmp_bi = where(tmp .gt. 365, tmp-365, tmp)
rclim1(:,:,:,:,1,:) = (/tmp_bi/)
tmp = rclim2(:,:,:,:,1,:)
tmp_tri = where(tmp .gt. 365, tmp-365, tmp)
rclim2(:,:,:,:,1,:) = (/tmp_tri/)
delete([/tmp,tmp_uni, tmp_bi, tmp_tri/])

;================================================================
filename = "Clim_CHIRPS.nc"
if (isfilepresent(filename)) then
	 system ("rm "+filename)
end if   
	f = addfile (filename, "c")

fAtt               = True            ; assign file attributes
   fAtt@title         = "East Africa Rainfall Onset Climatology"
   fAtt@source_file   = "climatology_chirps.ncl"
   fAtt@source_data   = "Filter_CHIRPS.nc"
   fAtt@Conventions   = "None"
   fAtt@creation_date = systemfunc ("date")
   fileattdef( f, fAtt )            ; copy file attributes

;=======DIMENSIONS==========================
filedimdef(f, "lon", dimsizes(lon), False)
filedimdef(f, "lat", dimsizes(lat), False)

seas_nr = ispan(0,2,1)
seas_nr!0 = "seas_nr"
seas_nr@name = "Season number"
filedimdef(f, "seas_nr", dimsizes(seas_nr), False)

convention = (/0, 1/)
convention!0 = "convention"
convention@name = "Convention relative(=0) and absolute(=1) pentad of the rainy season"
convention@comment = "Hydrological relative year starts at the driest month in the climatology = start_month"
filedimdef(f, "convention", 2, False)

pen = (/0, 1/)
pen!0 = "pen"
pen@name = "Pentad begin(=0) and end(=1) of the rainy season"
pen@long_name = "Rainy season in pentads since beginning of the hydrological year"
filedimdef(f, "pen", 2, False)

levels = (/5, 50, 95/)
levels!0 = "levels"
levels@name = "Confidence intervals"
filedimdef(f, "levels", 3, False)

;=========Variables=================================
threshold = new((/dimsizes(lat), dimsizes(lon), dimsizes(levels)/), "float", -9999.9)
threshold@name = "Rainfall threshold"
threshold!0 = "lat"
threshold&lat = lat
threshold!1 = "lon"
threshold&lon = lon
threshold!2 = "levels"
threshold&levels = levels
f->threshold = (/thres/)

ref_month = new((/dimsizes(lat), dimsizes(lon)/), "integer", -9999)
ref_month@name = "Month of lowest rainfall"
ref_month!0 = "lat"
ref_month&lat = lat
ref_month!1 = "lon"
ref_month&lon = lon
f->ref_month = (/ref_mon/)

clim = new(dimsizes(rclim), "integer", -9999)
clim@name = "Climatology of uni-modal rainfall"
clim!0 = "lat"
clim&lat = lat
clim!1 = "lon"
clim&lon = lon
clim!2 = "levels"
clim&levels = levels
clim!3 = "seas_nr"
clim&seas_nr = seas_nr
clim!4 = "convention"
clim&convention = convention
clim!5 = "pen"
clim&pen = pen
f->clim = (/rclim/)

clim1 = new(dimsizes(rclim), "integer", -9999)
clim1@name = "Climatology of bi-modal rainfall"
clim1!0 = "lat"
clim1&lat = lat
clim1!1 = "lon"
clim1&lon = lon
clim1!2 = "levels"
clim1&levels = levels
clim1!3 = "seas_nr"
clim1&seas_nr = seas_nr
clim1!4 = "convention"
clim1&convention = convention
clim1!5 = "pen"
clim1&pen = pen
f->clim1 = (/rclim1/)

clim2 = new(dimsizes(rclim), "integer", -9999)
clim2@name = "Climatology of tri-modal rainfall"
clim2!0 = "lat"
clim2&lat = lat
clim2!1 = "lon"
clim2&lon = lon
clim2!2 = "levels"
clim2&levels = levels
clim2!3 = "seas_nr"
clim2&seas_nr = seas_nr
clim2!4 = "convention"
clim&convention = convention
clim2!5 = "pen"
clim2&pen = pen
f->clim2 = (/rclim2/)

;; ============================================================

end
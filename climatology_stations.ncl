load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

; Load functions for some steps in this script
load "./functions.ncl"

begin
	file_in = "station_dual.nc"
  	in = addfile(file_in,"r")
	
	timeseries = in->precipitation_amount ; read rainfall data
	time = in->time
	
	ytmp = toint(cd_calendar(timeseries&time, 4)) ;years in timeseries of the file
	yy = ispan(1983,2013,1)
	dimYears = dimsizes(yy)-1

	ref_month = get_ref_mon(timeseries)	; Start month of hydrological year
  	;print(ref_month)
	
	t_opt=True
	t_opt@refmon = ref_month
	t_opt@nWeights = 9
	
	filtered_RR = new((/dimsizes(yy), 81/), "float", -999.9) 	; dimensions: number of years x (number of pentads(=72) + buffer for the filter(=9))
	do i=0,dimYears
		test = year_avail(timeseries, yy(i), t_opt) 	;get available data per year, buffer included
		avail = tofloat(num(.not. ismissing(timeseries({time |test(0):test(1)}))))/dimsizes(timeseries({time |test(0):test(1)}))
		if(avail .lt. 0.9)then		; 90% of the year should be available
			continue
		end if	
		if(.not.(ismissing(ind(time .eq. test(0))) .or. ismissing(ind(time .eq. test(1)))))then
			pentads = get_pentads(time, test, t_opt) 	;divide the year into pentads
			relrain = avg_rain(timeseries, pentads)		;get average rainfall per pentad
			filtered_RR(i,:) = get_filtering(relrain)		;filtering
		end if	
		delete(test)
		delete(avail)
	end do
	
	RR_new = dim_avg_n_Wrap(filtered_RR, 0) ; average pentad climatology
	nomiss = ind(.not. ismissing(RR_new))
	RR_perc = bootstrapping(filtered_RR(:,nomiss))	; get confidence intervals of pentad rainfall
	RR_1 = RR_new(nomiss)


bb = get_threshold(RR_1, True)	; get rainfall threshold
bb50 = bb[0]	; average threshold value

onset_clim50 = climatology_intervals(RR_1, bb50)	;get seasonality, onset and cessation pentads

bb5x = RR_perc(0,:)
bb95x = RR_perc(1,:)

;get seasonality, onset and cessation pentads for confidence intervals
bb5 = bb5x(bb[1])
onset_clim5 = climatology_intervals(RR_1, bb5)	

bb95 = bb95x(bb[1])
onset_clim95 = climatology_intervals(RR_1, bb95)

; convert onset/cessation dates to days since beginning of the year
ref_date = toint(cd_inv_calendar(1990,ref_month,1,0,0,0,"days since 1990-01-01",0))
onset_cal50 = 5*onset_clim50 + ref_date
onset_cal50 = where(onset_cal50 .gt. 365, onset_cal50-365, onset_cal50)

onset_cal5 = 5*onset_clim5 + ref_date
onset_cal5 = where(onset_cal5 .gt. 365, onset_cal5-365, onset_cal5)
onset_cal95 = 5*onset_clim95 + ref_date
onset_cal95 = where(onset_cal95 .gt. 365, onset_cal95-365, onset_cal95)

;=============================================================
; Write results to netcdf
filename = "climatology_file.nc"
if (isfilepresent(filename)) then
	system ("rm "+filename)
end if
   
f = addfile (filename, "c")

fAtt               = True            ; assign file attributes
  fAtt@title         = "Rainy Season Climatology"
  fAtt@source_file   = "climatology_stations.ncl"
  fAtt@source_data   = file_in
  fAtt@Conventions   = "None"
  fAtt@creation_date = systemfunc ("date")
  fileattdef( f, fAtt )            ; copy file attributes

;=======DIMENSIONS==========================
levels = (/5, 50, 95/)	
levels!0 = "levels"
filedimdef(f, "levels", 3, False)
	
	nq = max((/sum(dimsizes(onset_clim5)), sum(dimsizes(onset_clim50)), sum(dimsizes(onset_clim95))/)-2)
	print(nq)
	
	if(nq .gt. 1)then
		seas_nr = ispan(1,nq,1)
		seas_nr!0 = "seas_nr"
		seas_nr@name = "Season number"
		filedimdef(f, "seas_nr", nq, False)
	else	
		seas_nr = 1
		seas_nr!0 = "seas_nr"
		seas_nr@name = "Season number"
		filedimdef(f, "seas_nr", 1, False)
	end if	

convention = (/0, 1/)
convention!0 = "convention"
convention@name = "Convention relative(=0) and absolute(=1) pentad of the rainy season"
convention@comment = "Hydrological RELATIVE year starts at the driest month in the climatology = start_month"
filedimdef(f, "convention", 2, False)

pen = (/0, 1/)
pen!0 = "pen"
pen@name = "Begin(=0) and end(=1) of the rainy season"
pen@long_name = "Rainy season in pentads steps"
filedimdef(f, "pen", 2, False)

;=======VARIABLES==========================
threshold = new(dimsizes(levels), "float", -999.9) 
threshold!0 = "levels"
threshold&levels = levels
f->threshold = threshold

climatology = new((/dimsizes(levels), dimsizes(seas_nr), dimsizes(convention), dimsizes(pen)/), "integer", -999)
climatology!0 = "levels"
climatology&levels = levels
climatology!1 = "seas_nr"
climatology&seas_nr = seas_nr
climatology!2 = "convention"
climatology&convention = convention
climatology!3 = "pen"
climatology&pen = pen
climatology@start_month = ref_month
climatology@name = "Rainy season in pentad steps"
climatology@units = "Pentads = 5*days"
f->climatology = climatology

;-----------------------------
th  = f->threshold
if(all(ismissing(th)))then
	thres = (/bb5, bb50, bb95/)
	f->threshold = (/thres/)
	delete(th)
	delete(thres)
end if	

rs = f->climatology
seas_nr = rs&seas_nr


do i=0,dimsizes(seas_nr)-1
	tmp = all(ismissing(rs({levels| 5}, {seas_nr| seas_nr(i)}, {convention| :}, {pen| :})))
	print(tmp)
	if(dimsizes(dimsizes(onset_clim5)) .gt. 1)then
		if((i .lt. dimsizes(onset_clim5(:,0))) .and. (tmp .eq. True))then
			f->climatology({levels| 5}, {seas_nr| seas_nr(i)}, convention| 0, pen| :)=(/onset_clim5(i,:)/)
			f->climatology({levels| 5}, {seas_nr| seas_nr(i)}, convention| 1, pen| :)=(/onset_cal5(i,:)/)
		end if
	else
		if((i .lt. 1) .and. (tmp .eq. True))then
			f->climatology({levels| 5}, {seas_nr| seas_nr(i)}, convention| 0, pen| :)=(/onset_clim5/)
			f->climatology({levels| 5}, {seas_nr| seas_nr(i)}, convention| 1, pen| :)=(/onset_cal5/)
		end if
	end if	
	delete(tmp)
	
	tmp = all(ismissing(rs({levels| 50}, {seas_nr| seas_nr(i)}, {convention| :}, {pen| :})))
	print(tmp)
	if(dimsizes(dimsizes(onset_clim50)) .gt. 1)then
		if((i .lt. dimsizes(onset_clim50(:,0))) .and. (tmp .eq. True))then
			f->climatology({levels| 50}, {seas_nr| seas_nr(i)}, convention| 0, pen| :)=(/onset_clim50(i,:)/)
			f->climatology({levels| 50}, {seas_nr| seas_nr(i)}, convention| 1, pen| :)=(/onset_cal50(i,:)/)
		end if
	else
		if((i .lt. 1) .and. (tmp .eq. True))then
			f->climatology({levels| 50}, {seas_nr| seas_nr(i)}, convention| 0, pen| :)=(/onset_clim50/)
			f->climatology({levels| 50}, {seas_nr| seas_nr(i)}, convention| 1, pen| :)=(/onset_cal50/)
		end if
	end if	
	delete(tmp)
	
	tmp = all(ismissing(rs({levels| 95}, {seas_nr| seas_nr(i)}, {convention| :}, {pen| :})))
	print(tmp)
	if(dimsizes(dimsizes(onset_clim95)) .gt. 1)then
		if((i .lt. dimsizes(onset_clim95(:,0))) .and. (tmp .eq. True))then
			f->climatology({levels| 95}, {seas_nr| seas_nr(i)}, convention| 0, pen| :)=(/onset_clim95(i,:)/)
			f->climatology({levels| 95}, {seas_nr| seas_nr(i)}, convention| 1, pen| :)=(/onset_cal95(i,:)/)
		end if
	else
		if((i .lt. 1) .and. (tmp .eq. True))then
			f->climatology({levels| 95}, {seas_nr| seas_nr(i)}, convention| 0, pen| :)=(/onset_clim95/)
			f->climatology({levels| 95}, {seas_nr| seas_nr(i)}, convention| 1, pen| :)=(/onset_cal95/)
		end if
	end if	
	delete(tmp)
end do

end
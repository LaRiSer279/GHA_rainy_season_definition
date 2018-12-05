; Load functions for some steps in this script
load "./functions.ncl"

begin
; Read rainfall time series 
  	in = addfile("station_dual.nc","r")
	timeseries  = in->precipitation_amount ; Rainfall data
  
; Read climatology
	f = addfile ("climatology_file.nc", "r")
	
	;Get driest month 
	 mm1 = f->climatology@start_month 
	 print(mm1)
 	
 	t_opt=True
	t_opt@refmon = mm1
	t_opt@nWeights = 9
	ff = 1995
	
		;Check if enough data available
		
		; check start and end dates
 		test = year_avail(timeseries, ff, t_opt) 
 		; check the percentage of data
		avail = tofloat(num(.not. ismissing(timeseries({time |test(0):test(1)}))))/dimsizes(timeseries({time |test(0):test(1)}))
		
		if(avail .lt. 0.75)then
			print("period not available")
			exit
		end if	
		
		if(.not.(ismissing(ind(timeseries&time .eq. test(0))) .or. ismissing(ind(timeseries&time .eq. test(1)))))then
 			
 			;;;;; Split years into pentads ;;;;;;;;; 
			dd=get_pentads(timeseries&time, test, t_opt)
			;print(dd)
			
			;;;;;; Define rainfall onset ;;;;
			mmm = avg_rain(timeseries, dd) ; Calculate relative pentad rainfall
			mmm = where(ismissing(mmm) .eq. True, 0, mmm)
			
			
			mmmm = get_filtering(mmm)	; apply Lanczos Filter	
			print(mmmm)
			
			th  = f->threshold ; read climatological threshold values
			print(th)
			
			rs1 = f->climatology ; read climatological rainy seasons
			rs5 = rs1(0,:,0,:)	; rainy season(s) according to 5%-threshold
			rs50 = rs1(1,:,0,:)	; rainy season(s) according to average threshold
			rs95 = rs1(2,:,0,:)	; rainy season(s) according to 95%-threshold
			
			rp = get_rainy_periods(mmmm(ind(.not. ismissing(mmmm))), th(1))
			mmmmor50 = get_onset(rp,rs50(ind(.not. ismissing(rs50(:,0))),:)) 
			print(mmmmor50)
			
			rp5 = get_rainy_periods(mmmm(ind(.not. ismissing(mmmm))), th(0))
			mmmmor5 = get_onset(rp5,rs5(ind(.not. ismissing(rs5(:,0))),:)) 
			print(mmmmor5)
			
			rp95 = get_rainy_periods(mmmm(ind(.not. ismissing(mmmm))), th(2))
			mmmmor95 = get_onset(rp95,rs95(ind(.not. ismissing(rs95(:,0))),:)) 
			print(mmmmor95)
		end if	

end
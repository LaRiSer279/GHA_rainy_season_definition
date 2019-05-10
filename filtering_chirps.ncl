load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "functions.ncl"

begin
 file_in = "chirps_test.nc"
 in = addfile(file_in,"r")

;;;;;;============== Read in variables ================= ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  z  = in->precip ; Rainfall data
  lon = in->lon
  lat = in->lat

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ytmp = toint(cd_calendar(z&time, 4)) ;years in timeseries of the file
yrs = ispan(min(ytmp), max(ytmp),1)

;;;;;;;================ Pre-initialization =============;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ref_mon = new((/dimsizes(lat), dimsizes(lon)/), "integer", -999) ; SAVE
ref_mon@units = "months"
ref_mon!0 = "lat"
ref_mon&lat = lat
ref_mon!1 = "lon"
ref_mon&lon = lon

Tfilt_all = new((/dimsizes(lat), dimsizes(lon),dimsizes(yrs), 81/), "float", -999.9)
Tfilt_all!0 = "lat"
Tfilt_all&lat = lat
Tfilt_all!1 = "lon"
Tfilt_all&lon = lon
Tfilt_all!2 = "years"
Tfilt_all&years = yrs

Tfiltered = new((/dimsizes(lat), dimsizes(lon),dimsizes(yrs), 73/), "float", -999.9)
copy_VarCoords(Tfilt_all, Tfiltered)

;;;;;;;===================================================;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dimLon = dimsizes(lon)-1
dimLat = dimsizes(lat)-1  
dimYear = dimsizes(yrs)-1

do ll = 0,dimLon
	do la = 0,dimLat
  		timeseries = z(:,la,ll)
  
  		; skip grid points over water
  		if(all(ismissing(timeseries)))then
  			continue
  		end if	

		print("LON = " + lon(ll))
		print("LAT = " + lat(la))
	
		ref_mon(la,ll) = get_ref_mon(timeseries) ; get start of hydrological year
		print(ref_mon(la,ll))
			
		t_opt=True
		t_opt@refmon = ref_mon(la,ll)
		t_opt@nWeights = 9
	
		do i=0,dimYear
			; check availability
			test = year_avail(timeseries, yrs(i), t_opt)
			end_points = (/ind(timeseries&time .eq. test(0)), ind(timeseries&time .eq. test(1))/)
			
			if(all(.not. ismissing(end_points)) .eq. True)then
				subseries = timeseries({time |test(0):test(1)})
				avail = tofloat(num(.not. ismissing(subseries)))/dimsizes(subseries)
				
				if(avail .lt. 0.9)then
					delete(avail)
					delete(subseries)
					continue
				else
					; split timeseries into pentads
					pentads = get_pentads(subseries&time, test, t_opt)
					relrain = avg_rain(subseries, pentads)
					
					;filter
					Tfilt_all(la,ll,i,:) = get_filtering(relrain)
					
					; cut end points
					Tfiltered(la,ll,i,:) = Tfilt_all(la,ll,i,ind(.not. ismissing(Tfilt_all(la,ll,i,:))))
					delete(pentads)
					delete(relrain)
				end if
				delete(avail)
				delete(subseries)	
			end if	
			delete(test)
		end do
		delete(timeseries)
	end do
end do


;================================================================
filename = "Filter_CHIRPS.nc"
if (isfilepresent(filename)) then
	 system ("rm "+filename)
end if   
	f = addfile (filename, "c")

fAtt               = True            ; assign file attributes
   fAtt@title         = "East Africa Filtered Rainfall"
   fAtt@source_file   = "filtering_chirps.ncl"
   fAtt@source_data   =  file_in
   fAtt@Conventions   = "None"
   fAtt@creation_date = systemfunc ("date")
   fileattdef( f, fAtt )            ; copy file attributes

;=======DIMENSIONS==========================
filedimdef(f, "lon", dimsizes(lon), False)
filedimdef(f, "lat", dimsizes(lat), False)

years = yrs
years!0 = "years"
years@name = "years"
filedimdef(f, "years", dimsizes(years), False)

pentads = ispan(1,73,1)
pentads!0 = "pentads"
pentads@name = "pentads"
filedimdef(f, "pentads", dimsizes(pentads), False)

;=========Variables=================================
filtered_rainfall = new((/dimsizes(lat), dimsizes(lon), dimsizes(years), dimsizes(pentads)/), "float", -9999.9)
filtered_rainfall@name = "Filtered time series of rainfall"
filtered_rainfall!0 = "lat"
filtered_rainfall&lat = lat
filtered_rainfall!1 = "lon"
filtered_rainfall&lon = lon
filtered_rainfall!2 = "years"
filtered_rainfall&years = years
filtered_rainfall!3 = "pentads"
filtered_rainfall&pentads = pentads
f->filtered_rainfall = filtered_rainfall
f->filtered_rainfall = (/Tfiltered/)

ref_month = new((/dimsizes(lat), dimsizes(lon)/), "integer", -999)
ref_month@name = "Month of lowest rainfall"
ref_month!0 = "lat"
ref_month&lat = lat
ref_month!1 = "lon"
ref_month&lon = lon
f->ref_month = ref_month
f->ref_month = (/ref_mon/)			
end	
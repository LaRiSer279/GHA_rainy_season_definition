load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 


function get_ref_mon(timeseries[*]:numeric)

;	INPUT: rainfall timeseries including time as dimension
;	OUTPUT:	Month with the climatologically lowest average precipitation
;	CRITERION: min(running 3-month sum) 

local xtmp, ytmp, mtmp, ydim, y, m, ymprecip, minmon_val, tmp_min, d_tmp, tmp_ref
begin
	xtmp = cd_calendar(timeseries&time, 0) ;years in timeseries of the file
	ytmp = toint(xtmp(:,0))		;years
	mtmp = toint(xtmp(:,1))		;months
	ydim = tofloat(dimsizes(ispan(min(ytmp), max(ytmp), 1) ))	;number of available years
	ymprecip = new((/12/), "float", -999.9)	;array for testing
	
	do m = 0,11
		if(m .eq. 0)then
			ymprecip(m) = sum(timeseries(ind((mtmp .le. 2) .or. (mtmp .eq. 12))))
		else
			if(m .eq. 11)then
				ymprecip(m) = sum(timeseries(ind((mtmp .ge. 11) .or. (mtmp .eq. 1))))
			else		
				ymprecip(m) = sum(timeseries(ind((mtmp .ge. m) .and. (mtmp .le. m+2))))
			end if
		end if		
	end do
	;print(ymprecip)
	
	;index of minimum rainfall = month of minimum rainfall -1
	tmp_min = ind(ymprecip/ydim .eq. min(ymprecip/ydim))	
	if(dimsizes(tmp_min) .gt. 1)then	; case: multiple months detected
		; CRITERION: compute the sum of neighbour months and choose the lowest
		tmp_ref = new(dimsizes(tmp_min), "float", -999.9)
		d_tmp = dimsizes(tmp_min)-1
		do m=0,d_tmp
			if(tmp_min(m) .eq. 0)then
				tmp_ref(m) = ymprecip(11)+ymprecip(1)
			else
				if(tmp_min(m) .eq. 11)then
					tmp_ref(m) = ymprecip(10)+ymprecip(0)
				else
					tmp_ref(m) = ymprecip(tmp_min(m)-1)+ymprecip(tmp_min(m)+1)
				end if
			end if	
		end do
		minmon_val = tmp_min(ind(tmp_ref .eq. min(tmp_ref)))+1
	else							; case: one month detected
		minmon_val = tmp_min + 1
	end if	 
	
	; pass result to the main script
	return(minmon_val)
end
;-----------------------------------------------------------------------------------------
function year_avail(timeseries:numeric, year:integer, opt)
;	INPUT: 	- rainfall timeseries array including time as dimension
;			- year of interest
;			- options including: start month, number of weights for Lanczos-filtering
;	!!!!!!	currently only 9 weights are supported	!!!!!!!!
; 	
;	OUTPUT: array indices for the start and end points of the hydrological year of interest

local refmon, d_start, d_end, atts, buffer, offset
begin

atts = getvaratts(opt)	; get opt

;	get start month form opt
; 	if no start month is provided, determine start month
if((any(atts .eq. "refmon") .eq. True) .and. (ismissing(opt@refmon) .eq. False))then
	refmon = opt@refmon	 
else
	refmon = get_ref_mon_hist(timeseries)
end if

; determine buffer depending on time units and number of weights
if(isStrSubset(timeseries&time@units, "seconds"))then
	d_start = cd_inv_calendar(year,refmon,1,0,0,0,timeseries&time@units,0)
	d_start@units = timeseries&time@units
	;print(d_start@units)
	d_end = cd_inv_calendar(year+1,refmon,1,0,0,0,timeseries&time@units,0)
	d_end@units = timeseries&time@units
	
	if((any(atts .eq. "nWeights") .eq. True) .and. (ismissing(opt@nWeights) .eq. False))then
		buffer = 5*(opt@nWeights/2)
		;print(buffer)
		d_start = d_start - buffer*24*60*60
		;print(cd_calendar(d_start,-2))
		d_end = d_end + buffer*24*60*60-24*60*60
		;print(cd_calendar(d_end,-2))
	end if
else
	if(isStrSubset(timeseries&time@units, "days"))then
		d_start = cd_inv_calendar(year,refmon,1,0,0,0,timeseries&time@units,0)
		d_start@units = timeseries&time@units
		;print(d_start@units)
		d_end = cd_inv_calendar(year+1,refmon,1,0,0,0,timeseries&time@units,0)
		d_end@units = timeseries&time@units
		
		if((any(atts .eq. "nWeights") .eq. True) .and. (ismissing(opt@nWeights) .eq. False))then
			buffer = 5*(opt@nWeights/2)
			;print(buffer)
			d_start = d_start - buffer
			;print(cd_calendar(d_start,-2))
			d_end = d_end + buffer-1
			;print(cd_calendar(d_end,-2))
		end if	
	end if
end if		

; pass result to the main script
return((/d_start, d_end/))
end
;-----------------------------------------------------------------------------------------

function get_pentads(time[*]:numeric, period, opt)

; 	INPUT: dates, year, month of minimum rainfall, starting day of the year
; 	OUTPUT: indices for begin of pentades during the year
;	In case of leap years, Feb. 29th is added to 11th pentad, e.g. one year has ALWAYS 73 pentads

local tmp, pentads, p_start, p_end, cond, ld
begin
	
	p_start = ind(time .eq. period(0))	;starting point
	p_end = ind(time .eq. period(1))	; ending point
	tmp = ispan(p_start,p_end,5)		; span pentads
	
	; leap year correction
  	if(dimsizes(tmp) .gt. 72+opt@nWeights)then
  		cond = max(ind((isleapyear(cd_calendar(time(tmp), -1)/100) .eq. True) .and. (mod(cd_calendar(time(tmp), -1),100) .eq. 2)))
  		ld = where(tmp .gt. tmp(cond), tmp+1,tmp)
  		;print(ld)
  		pentads = ld(0:(dimsizes(tmp)-2))
  	else
  		pentads = tmp
  	end if
  	
    return(pentads)
end
;-----------------------------------------------------------------------------------------

function avg_rain(timeseries[*]:numeric, pentads[*]:integer)

;	INPUT: rainfall timeseries including time as dimension, Month of minimum rainfall, indices of pentads
;	OUTPUT: rainfall pentad averages 
	
local pmean, rank, p, tmp
begin
	rank = dimsizes(pentads)-1
	pmean = new(rank+1, "float", -999.9)

	do p=0,rank
		if(p .eq. rank)then
			if(any(timeseries(pentads(p):pentads(p)+4) .ne. timeseries@_FillValue) .eq. True)then
				tmp = timeseries(pentads(p):pentads(p)+4)
				pmean(p) = avg(tmp(ind(.not. ismissing(tmp))))
				delete(tmp)
			end if 	
		else
			if(any(timeseries(pentads(p):pentads(p+1)-1) .ne. timeseries@_FillValue) .eq. True)then
				tmp = timeseries(pentads(p):pentads(p+1)-1)
				pmean(p) = avg(tmp(ind(.not. ismissing(tmp))))
				delete(tmp)
			end if 
		end if
	end do
return(pmean)	
end 

;-----------------------------------------------------------------------------------------
function get_filtering(relTimeseries[*]:numeric)
; INPUT: rainfall pentad averages
; OUTPUT: filtered time series

local RR_new, sigma, ihp, nWgt, fca, wgt, nomissind, transf_relRR  
begin
	RR_new = new(dimsizes(relTimeseries), typeof(relTimeseries), -999.9)
	ihp = 0		; 	low-pass filtering
	sigma = 1	; 	power of the sigma factor, sigma = 1 is common
	
	nWgt = 9	; 	number of weights
	fca = 5/30.	;	cut-off frequency, in out case = 5 days per timestep/ 30 days 
	
	wgt = filwgts_lanczos(nWgt, ihp, fca, -999., sigma)		;	calculate filter weights
	RR_new  = wgt_runave(relTimeseries, wgt, 0)				; 	apply filter to time series
	
	return(RR_new)	; return results to main script 
end

;-----------------------------------------------------------------------------------------
function climatology_intervals(timeseries[*]:numeric, thres)
; INPUT: pentad time series, threshold value
; OUTPUT: climatological rainy seasons as start and end pentads

local onset_clim, rain_clim, emp, reduced, rind, i, length, o_tmp, dimclim, dim_clim, oclim
begin

rain_clim = ind(timeseries .ge. thres)	; indices when rainfall is above threshold
dim_clim = dimsizes(rain_clim)-1		

onset_clim = new((/10, 2/), "integer", -999)	; preset 10 rainy seasons
onset_clim(0,0) = rain_clim(0)
onset_clim(9,1) = rain_clim(dim_clim)

; get breaks between rainy periods
do i=1,dim_clim
	if(rain_clim(i)-rain_clim(i-1) .gt. 1)then
		emp = ind(ismissing(onset_clim(:,1)))
		onset_clim(emp(0),1) = rain_clim(i-1)
		emp = ind(ismissing(onset_clim(:,0)))
		onset_clim(emp(0),0) = rain_clim(i)
		delete(emp)
	end if
end do
;print(onset_clim)

; remove missing seasons
if(num(ismissing(onset_clim) .eq. True) .gt. 0)then
	oclim = dimsizes(onset_clim(:,0))-1
	do i=0,oclim
		reduced = ((ismissing(onset_clim(i,0)) .eq. False) .and. (ismissing(onset_clim(i,1)) .eq. True))
		if(reduced .eq. True)then
			onset_clim(i,1) = onset_clim(9,1)
			onset_clim(9,1) = onset_clim@_FillValue
			break
		end if
		delete(reduced)
	end do
	rind = onset_clim(ind(ismissing(onset_clim(:,0)) .eq. False),:)
	delete(onset_clim)
	onset_clim = rind
end if

; remove periods of less than 5 pentads
dimclim = dimsizes(onset_clim)
if(dimsizes(dimclim) .eq. 1)then
	o_tmp = onset_clim
else
	length = onset_clim(:,1)-onset_clim(:,0)
	if(num(length .ge. 5) .ge. 1)then
		o_tmp = onset_clim(ind(length .ge. 5),:)
		print(o_tmp)
	else
		o_tmp = onset_clim@_FillValue
	end if		
end if		

return(o_tmp)
end
;-----------------------------------------------------------------------------------------
function bootstrapping(timeseries[*][*]:numeric)

local nReps, rdim, tdim, i, j, dd1, bootind1, bootstrapStat1, bsdim, x5, x95, bb5x, bb95x, tdim_dim 
begin

nReps = 1000
rdim = nReps-1

tdim = dimsizes(timeseries)
tdim_dim = tdim(1)-1

bb5x = new(tdim(1), "float", -999.9)
bb95x = new(tdim(1), "float", -999.9)


	do j=0,tdim_dim
		dd1 = timeseries(:,j)
		bootstrapStat1 = new((/nReps/),"float")
		do i=0,rdim
			bootind1 = toint(ceil(random_uniform(0,1,tdim(0))*tdim(0))-1)
			;print(bootind1)
   			bootstrapStat1(i) = dim_avg(dd1(bootind1))
		end do
		qsort(bootstrapStat1)
		bsdim = dimsizes(bootstrapStat1)
		x5  = round(0.05*bsdim,3)-1
		x95 = round(0.95*bsdim,3)-1
	
		bb5x(j) = bootstrapStat1(x5)
		bb95x(j) = bootstrapStat1(x95)
		delete(dd1)
	end do

return((/bb5x, bb95x/))
end
;-----------------------------------------------------------------------------------------
function get_threshold(timeseries[*]:numeric, opt:logical)

local RR_1, tmp1, test2, thres, i, ttt, ttl, dtmp, ge8
begin

	RR_1 = timeseries(ind(.not. ismissing(timeseries)))
	tmp1 = dim_standardize(RR_1, 0)
	dtmp = dimsizes(tmp1)-2

	test2 = new((/2, dimsizes(tmp1)/), "float", -999)
	do i=1,dtmp
		test2(0,i) = i*((avg(tmp1(0:i)))^2)
		test2(1,i) = (dimsizes(tmp1)-1-i)*((avg(tmp1(i+1:dimsizes(tmp1)-1)))^2)
	end do
	test2(:,0) = 0
	test2(:,dimsizes(tmp1)-1) = 0
	ttt = local_max_1d((test2(1,:)+test2(0,:)), False, 0.0, 0)
	ttl = local_max_1d((test2(1,:)+test2(0,:)), False, 0.0, 1)
	;print_table([/RR_1,(test2(1,:)+test2(0,:))/], "%f, %f")
	;print(ttt)
	
	thres = [/ new(1, "float", -999.9), new(1, "integer", -999)/] 
	if(num(ttt .ge. 8.85) .ge. 1)then
		ge8 = ttl(ind(ttt .ge. 8.85))
		;print(RR_1(ge8))
		thres[0] = min(RR_1(ge8))
	else
		thres[0] = RR_1(ttl(ind(ttt .eq. max(ttt))))	
	end if
	
	if(opt .eq. True)then
		thres[1] = ind(RR_1 .eq. thres[0])
	end if
		
	return(thres)
end

;-----------------------------------------------------------------------------------------

function get_overlap(clim[2]:integer, series[2]:integer)
local cspan, sspan, overlap, i, cdim
begin

	cspan = ispan(clim(0),clim(1),1)
	sspan = ispan(series(0),series(1),1)
	
	cdim = dimsizes(sspan)-1
	overlap = 0
	do i=0,cdim
		if(any(sspan(i) .eq. cspan))then
			overlap = overlap + 1
		end if	
	end do

	return(overlap)
end 
;-----------------------------------------------------------------------------------------
function get_rainy_periods(trans_relRR[*]:numeric, threshold:float)
local rs_tmp, tmp, dim_rain, i, j, rperiods
begin

	rs_tmp = new((/30, 2/), "integer", -999)
	tmp = ind(trans_relRR .ge. threshold)
	
	dim_rain = dimsizes(tmp)-1
	
	j = 0
	do i=0,dim_rain
		if(dim_rain .eq. 0)then
			rs_tmp(j,0) = tmp(i)
			rs_tmp(j,1) = tmp(i)
			break			
		else
			if(dim_rain .eq. 1)then
				rs_tmp(j,0) = tmp(i)
				rs_tmp(j,1) = tmp(i+1)
				break
			else
				if(i .eq. 0)then
					rs_tmp(j,0) = tmp(i)
				else
					if(i .eq. dim_rain)then
						if(tmp(i) .gt. tmp(i-1)+1)then
							rs_tmp(j,1) = tmp(i-1)	
							j = j+1
							rs_tmp(j,0) = tmp(i)
							rs_tmp(j,1) = tmp(i)
						else	
							rs_tmp(j,1) = tmp(i)
						end if	
					else
						if(tmp(i) .gt. tmp(i-1)+1)then
							rs_tmp(j,1) = tmp(i-1)	
							j = j+1
							rs_tmp(j,0) = tmp(i)
						end if	
					end if
				end if
			end if	
		end if	
	end do
	

if(any( .not. ismissing(rs_tmp)) .eq. True)then
	rperiods = rs_tmp(ind(ismissing(rs_tmp(:,0)) .eq. False), :)
else
	rperiods = -999	
	rperiods@_FillValue = -999
end if

return(rperiods)

end
;-----------------------------------------------------------------------------------------
function get_onset(periods:numeric, onset_clim:integer)
local	p_dim, o_dim, isect, pent_isect, d_isect, i, j,isc, isc_tmp, rainy_season ,xx 
begin

if((all(ismissing(periods)) .eq. True) .or. (all(ismissing(onset_clim)) .eq. True))then
	print("one or more input is missing")
	exit
end if

p_dim = dimsizes(periods)
o_dim = dimsizes(onset_clim)
rainy_season = new(o_dim, "integer", -999)

;print(periods)
;print(onset_clim)
;print(p_dim)
if(dimsizes(p_dim) .gt. 1)then							;multiple pentad sequences
	if(dimsizes(o_dim) .eq. 1)then								;one rainy season
		;overlapping pentads per rainy season
		pent_isect = new(p_dim(0), "integer", -999)
		d_isect = dimsizes(pent_isect)-1
		 
		do i=0,d_isect
			pent_isect(i) = get_overlap(onset_clim, periods(i,:))
		end do
		
		if(any(pent_isect .gt. 0) .eq. True)then
			rainy_season(0) = min(periods(ind(pent_isect .gt. 0),0))
			rainy_season(1) = max(periods(ind(pent_isect .gt. 0),1))
		end if
	else												;multiple rainy seasons
		;overlapping pentads per rainy season
		pent_isect = new((/p_dim(0), o_dim(0)/), "integer", -999) 
		d_isect = dimsizes(pent_isect)-1
		
		;rainy season assignment
		isect = new(p_dim(0), "integer", -999) 

	
		do i=0,d_isect(0)
			do j=0,d_isect(1)
				pent_isect(i,j) = get_overlap(onset_clim(j,:), periods(i,:))
			end do
			if(any(pent_isect(i,:) .gt. 0) .eq. True)then
				;print(pent_isect(i,:))
				if(num(pent_isect(i,:) .eq. max(pent_isect(i,:))) .eq. 1)then
					isect(i) = ind(pent_isect(i,:) .eq. max(pent_isect(i,:)))
				else
					xx = ind(pent_isect(i,:) .eq. max(pent_isect(i,:)))
					isect(i) = xx(0)
				end if		
			end if	
		end do	
		;print(isect)
		
		do i=0,o_dim(0)-1
			if(any(.not. ismissing(ind(isect .eq. i))) .eq. True)then
				rainy_season(i,0) = min(periods(ind(isect .eq. i),0))
				rainy_season(i,1) = max(periods(ind(isect .eq. i),1))
			end if	
		end do
	end if
else													;one wet pentad sequence
	if(dimsizes(o_dim) .eq. 1)then								;one rainy season
		pent_isect = get_overlap(onset_clim, periods)
		if(pent_isect .gt. 0)then
			rainy_season(0) = min(periods)
			rainy_season(1) = max(periods)
		end if	
	else												;multiple rainy seasons
		pent_isect = new((/o_dim(0)/), "integer", -999)
		d_isect = dimsizes(pent_isect)-1
		
		do i=0,d_isect
			pent_isect(i) = get_overlap(onset_clim(i,:), periods)	;test for intersection with climatology
		end do
		;print(pent_isect)
		isect = ind(pent_isect .eq. max(pent_isect))
		;print(isect)
		if(any(.not. ismissing(isect)) .eq. True)then	;intersection 
			if(dimsizes(isect) .eq. 1)then
				rainy_season(isect,0) = min(periods)
				rainy_season(isect,1) = max(periods)
			else
				rainy_season(min(isect),0) = min(periods)
				rainy_season(min(isect),1) = min(periods)
			end if				
		end if	
	end if	
end if	

return(rainy_season)	
end
;-----------------------------------------------------------------------------------------

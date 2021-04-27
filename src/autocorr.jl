export sc_all, process_sc!, stream_autocorr

function sc_all(
	mseedfiles::AbstractArray,
	fs::AbstractFloat,
	cc_len::Real,
	cc_step::Real,
	freqmin::AbstractFloat,
	freqmax::AbstractFloat,
	maxlag::AbstractFloat,
	XMLDIR::String;
	responsefreq::Real=0.1,
)

	# read instrument response only once 
	RESP = read_resp(mseedfiles[1],XMLDIR)
	ind = findfirst(contains.(RESP.id,".BH"))
	RESP = RESP[ind]

	p = replace(basename(mseedfiles[1]),".ms"=>"") 
	pstr = join(split(p,"_",keepempty=false)," ")
	println("Auto-correlation $pstr $(now())")
	S = read_data("mseed",mseedfiles)

	process_sc!(
		S,
		fs,
		RESP,
		responsefreq=responsefreq,
	)
	F1 = seischannel2fft(S[1],cc_len,cc_step,freqmin,freqmax)
	F2 = seischannel2fft(S[2],cc_len,cc_step,freqmin,freqmax)
	F3 = seischannel2fft(S[3],cc_len,cc_step,freqmin,freqmax)

	C1 = correlate(F1,F2,maxlag)
	robuststack!(C1)
	C2 = correlate(F1,F3,maxlag)
	robuststack!(C2)
	C3 = correlate(F2,F3,maxlag)
	robuststack!(C3)
    return C1, C2, C3
end

function stream_autocorr(
	mseedfiles::AbstractArray,
	fs::AbstractFloat,
	cc_len::Real,
	cc_step::Real,
	freqmin::AbstractFloat,
	freqmax::AbstractFloat,
	maxlag::AbstractFloat,
	RESP::SeisData;
	responsefreq::Real=0.1,
	maxgaps::Int=100,
	aws_config::AWSConfig=global_aws_config(),
	bucket::String="scedc-pds",
	filetype::String="mseed",
	resp::Bool=true,
	t_max::Real=20.,
	factor::Real=10,
)
	# print net station channel and date 
	p = replace(mseedfiles[1],".ms"=>"") 
	println("Auto-correlation $p $(now())")

	# stream from S3 
	S = SeisData()
	try 
		if filetype == "mseed"
			S += mseed_stream(aws_config, bucket, mseedfiles)
		elseif filetype == "seisio"
			S += seisio_stream(aws_config, bucket, mseedfiles)
		end
	catch e
		julday = replace(basename(dirname(mseedfiles[1])),"_"=>"")
		d = yyyyjjj2date(julday)
		C = nancorr(d,fs,maxlag)
		return C,C,C
	end

	# short circuit if lots of gaps in data
	ngaps = get_gaps(S)
	if maximum(ngaps) > maxgaps
		julday = replace(basename(dirname(mseedfiles[1])),"_"=>"")
		d = yyyyjjj2date(julday)
		C = nancorr(d,fs,maxlag)
		return C,C,C
	end

	starttimes = [SeisIO.starttime(S.t[ii],S.fs[ii]) for ii = 1:S.n]
	endtimes = [SeisIO.endtime(S.t[ii],S.fs[ii]) for ii = 1:S.n]
	s = minimum(starttimes)
	t = maximum(endtimes) 

	# sync instrument response time 
	R = sync_resp(RESP,s,t)
	if length(R) == 0
		R = RESP[1]
	else
		R = R[1]
	end

	# remove instrument response
	process_sc!(
		S,
		fs,
		R,
		responsefreq=responsefreq,
		resp=resp,
		t_max=t_max,
		factor=factor,
	)

	if all(isnan.(S.x[1])) || length(S) != 3 
		julday = replace(basename(dirname(mseedfiles[1])),"_"=>"")
		d = yyyyjjj2date(julday)
		C = nancorr(d,fs,maxlag)
		return C,C,C
	end

	# convert to FFT 
	F1 = seischannel2fft(S[1],cc_len,cc_step,freqmin,freqmax)
	F2 = seischannel2fft(S[2],cc_len,cc_step,freqmin,freqmax)
	F3 = seischannel2fft(S[3],cc_len,cc_step,freqmin,freqmax)

	# correlate
	C1 = correlate(F1,F2,maxlag)
	robuststack!(C1)
	C2 = correlate(F1,F3,maxlag)
	robuststack!(C2)
	C3 = correlate(F2,F3,maxlag)
	robuststack!(C3)
    return C1,C2,C3
end

function process_sc!(
	S::SeisData,
	fs::Real,
	RESP::SeisChannel;
	responsefreq::Real=0.1,
	resp::Bool=true,
	t_max::Real=20.,
	factor::Real=10.,
)

	# check sampling rate for corrupted files
	S.fs = round.(S.fs)

	# if can't merge, return NaNs 
	try 
		merge!(S)
		sort!(S)
		taper!(S,t_max=t_max)
		ungap!(S,m=false)
		demean!(S)
		detrend!(S) 
	catch e
		for ii = 1:S.n
			S.x[ii] .= NaN
		end
		return S
	end

	# find and replace nans with zeros 
	nans = Dict()
	for ii = 1:S.n
		nans[ii] = findall(isnan.(S.x[ii]))
		S.x[ii][nans[ii]] .= 0 
	end

	# highpass filter 
	highpass!(S,responsefreq,zerophase=true,corners=2)

	# replace gaps with zeros
	for ii = 1:S.n
		S.x[ii][nans[ii]] .= 0 
	end 

	# resample 
	resample!(S,fs=fs) 

	# find and replace nans with zeros 
	nans = Dict()
	for ii = 1:S.n
		nans[ii] = findall(isnan.(S.x[ii]))
		S.x[ii][nans[ii]] .= 0 
	end
      
	# phase shift if necessary
	allt = vcat(S.t...)[1:2:end,2]
	if all(t->t==allt[1],allt)
		ϕshift = false
	else
		ϕshift = true
	end
    phase_shift!(S, ϕshift=ϕshift)

	# remove instrument response
	if resp
		for ii = 1:S.n
			S.loc[ii] = RESP.loc
			S.gain[ii] = RESP.gain
			S.resp[ii] = RESP.resp
		end
		remove_resp!(S)
	end

	# setup taper windows 
	winsamples = round(Int,t_max * S.fs[1])
	winsamples += iseven(winsamples)
	window = SeisIO.hanning(2 * winsamples - 1)
	leftwin = window[winsamples:end]
	rightwin = window[1:winsamples]

	# find outlier points and taper 
	for ii = 1:S.n 
		# find points 10x greater than MAD 
		fm = fast_mad(S.x[ii], round(Int,length(S.x[ii]) / 100)) * factor 
		ind = findall(abs.(S.x[ii]) .> fm)
		if length(ind) <= 1 
			continue
		end

		# # get smooth version of envelope
		# eup, edown = envelope(S.x[ii])
		# smooth!(eup,round(Int,S.fs[1] / 2))

		# # get new fm based on envelope 
		# fm = fast_mad(eup) 
		# ind = findall(eup .> 2 * fm)

		# if length(ind) <= 1 
		# 	continue
		# end

		leftind = [ind[1]; ind[findall(diff(ind) .> 1) .+ 1]]
		rightind = [ind[findall(diff(ind) .> 1)]; ind[end]]
		
		# taper high amplitudes on the left 
		for jj = 1:length(rightind)
			if rightind[jj] + winsamples + 1 > length(S.x[ii])
				continue
			end
			S.x[ii][rightind[jj]+1:rightind[jj]+winsamples] .*= rightwin
		end

		# taper high amplitudes on the right 
		for jj = 1:length(leftind)
			if leftind[jj] - winsamples < 1
				continue
			end
			S.x[ii][leftind[jj]-winsamples:leftind[jj]-1] .*= leftwin
		end

		# set all points between right and left indices to zero 
		for jj = 1:length(leftind)
			S.x[ii][leftind[jj]:rightind[jj]] .= 0 
		end
	end

	# replace gaps with zeros
	for ii = 1:S.n
		S.x[ii][nans[ii]] .= 0 
	end

	return nothing
end

function seischannel2fft(
	C::SeisChannel,
	cc_len::Real,
	cc_step::Real,
	freqmin::Real,
	freqmax::Real;
)

	# convert to RawData
	R = RawData(C,cc_len,cc_step) 
	ind = findall(abs.(R.x) .< fast_mad2d(R.x) ./ 1000)
	detrend!(R)
	taper!(R)
	# highpass!(R,freqmin,zerophase=true)
	# whiten!(R,freqmin,freqmax)
	R.x[ind] .= 0
	onebit!(R)
	return rfft(R)
end
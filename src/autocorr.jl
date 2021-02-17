export sc_all, process_sc!, stream_autocorr

function sc_all(
	mseedfiles::AbstractArray,
	fs::AbstractFloat,
	cc_len::Real,
	cc_step::Real,
	freqmin::AbstractFloat,
	freqmax::AbstractFloat,
	maxlag::AbstractFloat,
	CORROUT::String,
	XMLDIR::String;
	responsefreq::Real=0.1,
)

	# read instrument response only once 
	RESP = read_resp(mseedfiles[1],XMLDIR)
	ind = findfirst(contains.(RESP.id,".HH"))
	RESP = RESP[ind]

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
	stack!(C1)
	save_corr(cpu(C1),CORROUT)
	C2 = correlate(F1,F3,maxlag)
	stack!(C2)
	save_corr(cpu(C2),CORROUT)
	C3 = correlate(F2,F3,maxlag)
	stack!(C3)
	save_corr(cpu(C3),CORROUT)
    return nothing
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
)
	# print net station channel and date 
	p = replace(basename(mseedfiles[1]),".ms"=>"") 
	pstr = join(split(p,"_",keepempty=false)," ")
	println("Auto-correlation $pstr $(now())")

	# stream from S3 
	S = SeisData()
	for ii = 1:3 
		S += SCEDC.s3_get_seed(
			"scedc-pds",
			mseedfiles[ii],
			false,
			false,
			false,
			false,
			false,
			false,
			false,
			false,
			false,
			0,
		)
	end

	# short circuit if lots of gaps in data
	ngaps = get_gaps(S)
	if maximum(ngaps) > maxgaps
		julday = replace(basename(dirname(mseedfiles[1])),"_"=>"")
		d = yyyyjjj2date(julday)
		C = nancorr(S,d,fs,maxlag)
		return C,C,C
	end

	starttimes = [SeisIO.starttime(S.t[ii],S.fs[ii]) for ii = 1:S.n]
	endtimes = [SeisIO.endtime(S.t[ii],S.fs[ii]) for ii = 1:S.n]
	s = minimum(starttimes)
	t = maximum(endtimes) 

	# sync instrument response time 
	R = sync_resp(RESP,s,t)
	if length(RESP) == 0
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
	)

	# convert to FFT 
	F1 = seischannel2fft(S[1],cc_len,cc_step,freqmin,freqmax)
	F2 = seischannel2fft(S[2],cc_len,cc_step,freqmin,freqmax)
	F3 = seischannel2fft(S[3],cc_len,cc_step,freqmin,freqmax)

	# correlate
	C1 = correlate(F1,F2,maxlag)
	stack!(C1)
	C2 = correlate(F1,F3,maxlag)
	stack!(C2)
	C3 = correlate(F2,F3,maxlag)
	stack!(C3)
    return C1,C2,C3
end

function process_sc!(
	S::SeisData,
	fs::Real,
	RESP::SeisChannel;
	responsefreq::Real=0.1,
)
	merge!(S)
	sort!(S)
	taper!(S,t_max=100.)
    ungap!(S,m=false)
	demean!(S)
	detrend!(S) 

	# find and replace nans with zeros 
	nans = Dict()
	for ii = 1:3
		nans[ii] = findall(isnan.(S.x[ii]))
		S.x[ii][nans[ii]] .= 0 
	end

	# highpass filter 
	highpass!(S,responsefreq,zerophase=true,corners=2)
      
	# phase shift if necessary
	if S[1].t[1,2] == S[2].t[1,2] == S[3].t[1,2]
		ϕshift = false
	else
		ϕshift = true
	end
    phase_shift!(S, ϕshift=ϕshift)

	# remove instrument response
	for ii = 1:3
		S.loc[ii] = RESP.loc
		S.gain[ii] = RESP.gain
		S.resp[ii] = RESP.resp
	end
	remove_resp!(S)
	resample!(S,fs=fs) 
	
	# replace gaps with zeros
	for ii = 1:3
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
	detrend!(R)
	taper!(R)
	bandpass!(R,freqmin,freqmax,zerophase=true)
	clip!(R,3)
	F = rfft(R)
	whiten!(F,freqmin,freqmax)
	return F
end
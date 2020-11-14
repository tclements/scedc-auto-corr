export sc_all, process

function sc_all(
	files::AbstractArray,
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
	println("Reading $(basename(dirname(files[1]))) $(replace(basename(files[1])[1:7],"_"=>""))")

	# read instrument response only once 
	RESP = read_resp("sxml",instpath,s=s,t=t)

	F1 = process_sc(
		files[1],
		fs,
		RESP,
		cc_len,
		cc_step,
		freqmin,
		freqmax,
		responsefreq=responsefreq,
	)

	# read 2nd channel
	F2 = process_sc(
		files[2],
		fs,
		RESP,
		cc_len,
		cc_step,
		freqmin,
		freqmax,
		responsefreq=responsefreq,
	)

	# read 3rd channel
	F3 = process_sc(
		files[3],
		fs,
		RESP,
		cc_len,
		cc_step,
		freqmin,
		freqmax,
		responsefreq=responsefreq,
	)

	C1 = correlate(F1,F2,maxlag)
	stack!(C1)
	save_corr(C1,CORROUT)
	C2 = correlate(F1,F3,maxlag)
	stack!(C2)
	save_corr(C2,CORROUT)
	C3 = correlate(F2,F3,maxlag)
	stack!(C3)
	save_corr(C3,CORROUT)
    return nothing
end

function process_sc(
	file::String,
	fs::Real,
	RESP::SeisData,
	cc_len::Real,
	cc_step::Real,
	freqmin::Real,
	freqmax::Real;
	responsefreq::Real=0.1,
)
	S = read_data("mseed",file)
	merge!(S)
    ungap!(S)
	detrend!(S)         # remove mean & trend from channel
	taper!(S)                      # taper channel ends
	resample!(S,fs=fs) # downsample to lower fs
	taper!(S)
    phase_shift!(S, ϕshift=true) # t

	# remove instrument response
	net,sta,loc,chan = split(S.id[1],'.')
	ind = findfirst(RESP.id .== S[1].id)
	highpass!(S,responsefreq,zerophase=true,corners=2)
	S.loc[1] = RESP[ind].loc
	S.gain[1] = RESP[ind].gain
	S.resp[1] = RESP[ind].resp
	remove_resp!(S)

	# convert to RawData
	R = RawData(S,cc_len,cc_step)
	detrend!(R)
	taper!(R)
	bandpass!(R,freqmin,freqmax,zerophase=true)
	clip!(R,3)
	F = rfft(R)
	whiten!(F,freqmin,freqmax)
	return F
end
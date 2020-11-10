module SCEDCAutoCorr


using SeisIO,SeisNoise, AWSS3, AWSCore, Dates, Glob, Statistics
export sc_all, prunefiles,upload_par, process

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
	s = date_yyyyddd(files[1][end-9:end-3])
	t = s + Day(1)
	s = Dates.format(s, "yyyy-mm-dd HH:MM:SS")
	t = Dates.format(t, "yyyy-mm-dd HH:MM:SS")
	net = basename(files[1])[1:2]
	sta = split(basename(files[1]),"_")[1][3:end]
	instpath = joinpath(XMLDIR,net * '_' * sta * ".xml" )
	RESP = read_meta("sxml",instpath,s=s,t=t)

	F1 = process(
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
	F2 = process(
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
	F3 = process(
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

function process(
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

function prunefiles(filelist::AbstractArray)
    if length(filelist) == 0
	    return []
    end
	files = deepcopy(filelist)
    fsizes = [filesize(f) for f in files]
    ind = findall((fsizes .< median(fsizes) / 4) .| (fsizes .> median(fsizes) * 2))
	deleteat!(files,ind)

    # get individual stations
    stations = unique([replace(basename(f)[1:7],"_"=>"") for f in files])
    infiles = []
    for jj in eachindex(stations)
	    ind = findall(occursin.(stations[jj],files))
        if length(ind) != 3
	        continue
        end
        push!(infiles,files[ind])
    end
	return infiles
end

function date_yyyyddd(yearday::String)
    @assert occursin(r"[1-2][0-9][0-9][0-9][0-3][0-6][0-9]",yearday)
    yint = parse(Int,yearday[1:4])
    dint = parse(Int,yearday[5:end])
    @assert dint <= 366 "Input day must be less than or equal to 366"
    return DateTime(yint) + Day(dint-1)
end

function upload_par(aws::Dict,output_bucket::String,s3file::String,ec2file::String)
    println("Uploading file $ec2file")
    s3_put(aws,output_bucket,s3file,read(ec2file))
end
end

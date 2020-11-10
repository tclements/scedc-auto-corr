module SCEDCAutoCorr


using SeisIO,SeisNoise, AWSS3, AWSCore, Dates, Glob, Statistics
export sc_all, prunefiles,upload_par, process

function sc_all(files::AbstractArray,fs::AbstractFloat,cc_len::Int,cc_step::Int,
	            freqmin::AbstractFloat,freqmax::AbstractFloat,maxlag::AbstractFloat,
				CORROUT::String,XMLDIR::String;maxgaps::Int=10)
	println("Reading $(basename(dirname(files[1]))) $(replace(basename(files[1])[1:7],"_"=>""))")
	try
		F1 = process(files[1],fs,XMLDIR, cc_len, cc_step, freqmin, freqmax,maxgaps=maxgaps)

		# read 2nd channel
		F2 = process(files[2],fs,XMLDIR, cc_len, cc_step, freqmin, freqmax,maxgaps=maxgaps)

		# read 3rd channel
		F3 = process(files[3],fs,XMLDIR, cc_len, cc_step, freqmin, freqmax,maxgaps=maxgaps)

		C1 = correlate(F1,F2,maxlag)
		stack!(C1)
		save_corr(C1,CORROUT)
		C2 = correlate(F1,F3,maxlag)
		stack!(C2)
		save_corr(C2,CORROUT)
		C3 = correlate(F2,F3,maxlag)
		stack!(C3)
		save_corr(C3,CORROUT)
	catch e 
		println(e)
	end
    return nothing
end

function process(
	file::String,
	fs::Real,
	XMLDIR::String,
	cc_len::Real,
	cc_step::Real,
	freqmin::Real,
	freqmax::Real;
	maxgaps::Int=10,
)
	S = read_data("mseed",file)
	merge!(S)
    ungap!(S)
	detrend!(S)         # remove mean & trend from channel
	taper!(S)                      # taper channel ends
	if fs ∉ S.fs
		filtfilt!(S,fh=fs/2,rt="Lowpass")    # lowpass filter before downsampling
	end
	resample!(S,fs=fs) # downsample to lower fs
	taper!(S)
    phase_shift!(S, ϕshift=true) # t

	# remove instrument response
	start = Dates.unix2datetime(S[1].t[1,2] * 1e-6)
	s = string(round(start,Dates.Day))
	t = string(round(start,Dates.Day) + Day(1))
	net,sta,loc,chan = split(S.id[1],'.')
	instpath = joinpath(XMLDIR,net * '_' * sta * ".xml" )
	RESP = read_meta("sxml",instpath,s=s,t=t)
	ind = findfirst(RESP.id .== S[1].id)
	highpass!(S,0.1,zerophase=true,corners=2)
	S.loc[1] = RESP[ind].loc
	S.gain[1] = RESP[ind].gain
	translate_resp!(S,RESP[ind].resp)
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
    ind = findall(fsizes .< median(fsizes) / 10)
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

function upload_par(aws::Dict,output_bucket::String,s3file::String,ec2file::String)
    println("Uploading file $ec2file")
        s3_put(aws,output_bucket,s3file,read(ec2file))
	end
end

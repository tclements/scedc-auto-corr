__precompile__()
module AC

using SeisIO,SeisNoise, JLD2, AWSS3, AWSCore, Dates, Glob, Statistics
export sc_all, scedc_files,upload_par

function sc_all(files::AbstractArray,fs::AbstractFloat,cc_len::Int,cc_step::Int,
	            freqmin::AbstractFloat,freqmax::AbstractFloat,maxlag::AbstractFloat,
				CORROUT::String;maxgaps::Int=10)
	println("Reading $(basename(dirname(files[1]))) $(replace(basename(files[1])[1:7],"_"=>""))")
	try
		S1 = read_data("mseed",files[1])
		if size(S1[1].t,1) > maxgaps
			return nothing
		end
		process_raw!(S1,fs)

		# remove instrument response
 	    start = Dates.unix2datetime(S1[1].t[1,2] * 1e-6)
	    s = string(round(start,Dates.Day))
		t = string(round(start,Dates.Day) + Day(1))
		net,sta,loc,chan = split(S1.id[1],'.')
		instpath = joinpath(expanduser("~/XML"),string(year(start)),net * '_' * sta * ".xml" )
		RESP = read_sxml(instpath,s=s,t=t)
		ind = findfirst(RESP.id .== S1[1].id)
		highpass!(S1,0.01,zerophase=true)
		S1.loc[1] = RESP[ind].loc
		S1.gain[1] = RESP[ind].gain
		translate_resp!(S1,RESP[ind].resp)

		# convert to RawData
	    R1 = RawData(S1,cc_len,cc_step)
		S1 = nothing
		detrend!(R1)
		taper!(R1)
		bandpass!(R1,freqmin,freqmax,zerophase=true)
		clip!(R1,3)
		F1 = compute_fft(R1)
		R1 = nothing
		whiten!(F1,freqmin,freqmax)

		# read 2nd channel
		S2 = read_data("mseed",files[2])
		if size(S2[1].t,1) > maxgaps
		return nothing
		end
		process_raw!(S2,fs)
		ind = findfirst(RESP.id .== S2[1].id)
		highpass!(S2,0.01,zerophase=true)
		S2.loc[1] = RESP[ind].loc
		S2.gain[1] = RESP[ind].gain
		translate_resp!(S2,RESP[ind].resp)
		R2 = RawData(S2,cc_len,cc_step)
		S2 = nothing
		detrend!(R2)
		taper!(R2)
		bandpass!(R2,freqmin,freqmax,zerophase=true)
		clip!(R2,3)
		F2 = compute_fft(R2)
		R2 = nothing
		whiten!(F2,freqmin,freqmax)

		# read 3rd channel
		S3 = read_data("mseed",files[3])
		if size(S3[1].t,1) > maxgaps
			return nothing
		end
		process_raw!(S3,fs)
		ind = findfirst(RESP.id .== S3[1].id)
		highpass!(S3,0.01,zerophase=true)
		S3.loc[1] = RESP[ind].loc
		S3.gain[1] = RESP[ind].gain
		translate_resp!(S3,RESP[ind].resp)
		R3 = RawData(S3,cc_len,cc_step)
		S3 = nothing
		detrend!(R3)
		taper!(R3)
		bandpass!(R3,freqmin,freqmax,zerophase=true)
		clip!(R3,3)
		F3 = compute_fft(R3)
		R3 = nothing
		whiten!(F3,freqmin,freqmax)

		C1 = compute_cc(F1,F2,maxlag)
		if isnothing(C1)
			return nothing
		end
		stack!(C1,stacktype=robuststack)
		save_corr(C1,CORROUT)
		C2 = compute_cc(F1,F3,maxlag)
		if isnothing(C2)
			return nothing
		end
		stack!(C2,stacktype=robuststack)
		save_corr(C2,CORROUT)
		C3 = compute_cc(F2,F3,maxlag)
		if isnothing(C3)
			return nothing
		end
		stack!(C3,stacktype=robuststack)
		save_corr(C3,CORROUT)
	catch
	end
    return nothing
end

function scedc_files(DATA::String,date::Date)
	filedir = joinpath(DATA,SeisNoise.scedcpath(date))
    files = glob("*",filedir)
    if length(files) == 0
	    return []
    end
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

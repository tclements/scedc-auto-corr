export prunefiles, upload_par, read_resp, date_yyyyddd

function prunefiles(filelist::AbstractArray; minsize = 0.25, maxsize = 2.)
    if length(filelist) == 0
	    return []
    end
	files = deepcopy(filelist)
    fsizes = [filesize(f) for f in files]
    ind = findall((fsizes .< median(fsizes) * minsize) .| (fsizes .> median(fsizes) * maxsize))
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

function read_resp(file::String,XMLDIR::String)
    s = date_yyyyddd(file[end-9:end-3])
	t = s + Day(1)
	s = Dates.format(s, "yyyy-mm-dd HH:MM:SS")
	t = Dates.format(t, "yyyy-mm-dd HH:MM:SS")
	net = basename(file)[1:2]
	sta = split(basename(file),"_")[1][3:end]
	instpath = joinpath(XMLDIR,net * '_' * sta * ".xml" )
    return read_meta("sxml",instpath,s=s,t=t)
end
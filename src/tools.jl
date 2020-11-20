export prunefiles, upload_par, read_resp, yyyyjjj2date, date2yyyyjjj, XML_download, indexpath, scedcpath

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

function yyyyjjj2date(yearday::String)
    @assert occursin(r"[1-2][0-9][0-9][0-9][0-3][0-9][0-9]",yearday)
    yint = parse(Int,yearday[1:4])
    dint = parse(Int,yearday[5:end])
    @assert dint <= 366 "Input day must be less than or equal to 366"
    return DateTime(yint) + Day(dint-1)
end

function date2yyyyjjj(d::TimeType)
    return "$(year(d))_$(lpad(dayofyear(d),3,"0"))"
end

function upload_par(aws::Dict,output_bucket::String,s3file::String,ec2file::String)
    println("Uploading file $ec2file")
    s3_put(aws,output_bucket,s3file,read(ec2file))
end

function read_resp(file::String,XMLDIR::String)
    s = yyyyjjj2date(file[end-9:end-3])
	t = s + Day(1)
	s = Dates.format(s, "yyyy-mm-dd HH:MM:SS")
	t = Dates.format(t, "yyyy-mm-dd HH:MM:SS")
	net = basename(file)[1:2]
	sta = split(basename(file),"_")[1][3:end]
	instpath = joinpath(XMLDIR,net * '_' * sta * ".xml" )
    return read_meta("sxml",instpath,s=s,t=t)
end

function XML_download(aws,XMLDIR)
    if !isdir(XMLDIR)
        mkpath(XMLDIR)
    end
    req = collect(s3_list_objects(aws,"scedc-pds","FDSNstationXML/CI/"))
    xmlin = [r["Key"] for r in req]
    xmlout = joinpath.(XMLDIR,basename.(xmlin))
    for ii = 1:length(xmlin)
        s3_get_file(aws,"scedc-pds",xmlin[ii],xmlout[ii])
    end
    return nothing
end

function indexpath(d::Date)
    days = (d - Date(Year(d))).value + 1
    n = ndigits(days)
	jstr = ('0' ^ (3 - n)) * string(days)
	ystr = string(Year(d).value)
    outstring = "continuous_waveforms/index/csv/year="
    outstring *= ystr * "/year_doy="
    outstring *= ystr * '_' * jstr
    outstring *= "/$(ystr)_$(jstr)_waveform_index.csv"
    return outstring
end

"""
  scedcpath(filename)
Convert filename to scedc-pds path.
"""
function scedcpath(filename::String)
    year = filename[14:17]
    day = filename[18:20]
    return "continuous_waveforms/" * year * '/' * year * '_' * day * '/' * filename
end
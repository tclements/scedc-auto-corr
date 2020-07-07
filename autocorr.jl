YYYY = parse(Int,ARGS[1])
using Distributed
addprocs()
@eval @everywhere YYYY=$YYYY
@everywhere push!(LOAD_PATH,pwd())
@everywhere using AC
@everywhere begin
using SeisIO, SeisNoise,Glob, Dates, AWSCore, AWSS3

# file cleaning
startdate = Date(YYYY,1,1)
enddate = Date(YYYY+1,1,1) - Day(1)
channel = "BH*"
network = "CI"
input_bucket = "scedc-pds"
output_bucket = "scedc-auto-corr"
aws = aws_config(region = "us-west-2")
date_range = startdate:Day(1):enddate
fs = 20.
cc_len, cc_step = 1800, 450
freqmin = 0.1
freqmax = 9.5
maxlag = 60.

# create dirs
DATA = joinpath(expanduser("~/data/"),string(YYYY))
CORR = joinpath(expanduser("~/CORR"),string(YYYY))
XML = joinpath(expanduser("~/XML"),string(YYYY))

if !isdir(DATA)
    mkpath(DATA)
end
if !isdir(CORR)
    mkpath(CORR)
end

if !isdir(XML)
    mkpath(XML)
end

# get stationXML
XMLlist = s3_get(aws,output_bucket,"XMLlist.txt")
xmlfiles = split(XMLlist,"\n")[1:end-1]
s3XML = [joinpath("FDSNstationXML/CI",x) for x in xmlfiles]
ec2XML = [joinpath(XML,x) for x in xmlfiles]
Nxml = length(s3XML)
end

pmap(s3_get_file,fill(aws,Nxml),fill(input_bucket,Nxml),s3XML,ec2XML)
t1 = now()

for ii in eachindex(date_range)
    date = date_range[ii]
    scedctransfer(DATA,date,channel=channel,network=network)
    infiles = scedc_files(DATA,date)
    N = length(infiles)
    if length(infiles) == 0
        continue
    end

    pmap(sc_all,infiles,fill(fs,N),fill(cc_len,N),fill(cc_step,N),fill(freqmin,N),
    fill(freqmax,N),fill(maxlag,N),fill(CORR,N))
    rm.(collect(Iterators.flatten(infiles)))
end
t2 = now()
println("Computation took $(t2-t1)")

# upload files to s3
ec2files = glob("*.jld2",CORR)
s3files = [replace(s,"/home/ubuntu/"=>"") for s in ec2files]
Ncorr = length(ec2files)
@eval @everywhere s3files=$s3files
@eval @everywhere ec2files=$ec2files
pmap(upload_par,fill(aws,Ncorr),fill(output_bucket,Ncorr),s3files,ec2files)

module SCEDCCorr

using SeisIO, SeisNoise, AWSS3, AWSCore, Dates, Glob, Statistics, JLD2
using SCEDC, CSV, DataFrames
include("tools.jl")
include("autocorr.jl")
include("LHcorrelation.jl")
end
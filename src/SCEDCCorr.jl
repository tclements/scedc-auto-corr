module SCEDCCorr

using SeisIO, SeisNoise, AWSS3, AWS, Dates, Glob, Statistics, JLD2
using CSV, DataFrames, SCEDC, Serialization
include("tools.jl")
include("autocorr.jl")
include("LHcorrelation.jl")
end
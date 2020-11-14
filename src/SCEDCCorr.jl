module SCEDCCorr

using SeisIO,SeisNoise, AWSS3, AWSCore, Dates, Glob, Statistics
include("tools.jl")
include("autocorr.jl")
include("LHcorrelation.jl")
end
export prepare_LH, LH_to_FFT, read_and_remove, all2all!

function prepare_LH(
    fileE::String,
    fileN::String,
    fileZ::String,
    XMLDIR::String,
    freqmin::Real,
    freqmax::Real,
    cc_len::Real,
    cc_step::Real,
)
    RESP = read_resp(fileE,XMLDIR)

    # read and remove instrument resp
    E = read_and_remove(fileE,freqmin,freqmax,RESP)
    N = read_and_remove(fileN,freqmin,freqmax,RESP)
    Z = read_and_remove(fileZ,freqmin,freqmax,RESP)

    # bandpass between 15 and 50 seconds 
    Efilt = filtfilt(E,rt="Bandpass",fl=1/50.,fh=1/15.,np=2)
    Nfilt = filtfilt(E,rt="Bandpass",fl=1/50.,fh=1/15.,np=2)
    Zfilt = filtfilt(E,rt="Bandpass",fl=1/50.,fh=1/15.,np=2)

    # normalize Z comp with 128-second weight 
    weights = smooth(abs.(Zfilt.x[1]),64)
    Z.x[1] ./= weights 

    # normalize E and N comps with 128-second weight 
    Eweights = smooth(abs.(Efilt.x[1]),64)
    Nweights = smooth(abs.(Nfilt.x[1]),64)
    weights = maximum(hcat([Eweights,Nweights]...),dims=2)[:]
    E.x[1] ./= weights
    N.x[1] ./= weights 

    # convert to FFT 
    return LH_to_FFT(E,N,Z,cc_len,cc_step)
end

function LH_to_FFT(E,N,Z,cc_len,cc_step)
    Eraw = RawData(E,cc_len,cc_step)
    Nraw = RawData(N,cc_len,cc_step)
    Zraw = RawData(Z,cc_len,cc_step)

    # simple processing 
    detrend!(Eraw)
    detrend!(Nraw)
    detrend!(Zraw)
    taper!(Eraw,max_length=100.)
    taper!(Nraw,max_length=100.)
    taper!(Zraw,max_length=100.)

    # convert to fft 
    Efft = rfft(Eraw)
    Nfft = rfft(Nraw)
    Zfft = rfft(Zraw)

    # whiten E,N,Z components with 100 s amplitude smoothing 
    Efft.fft ./= smooth(abs.(Efft.fft),50)
    Nfft.fft ./= smooth(abs.(Efft.fft),50)
    Zfft.fft ./= smooth(abs.(Zfft.fft),50)
    return Efft, Nfft, Zfft 
end

function read_and_remove(file::String,freqmin::Real,freqmax::Real,RESP::SeisData)
    S = read_data("mseed",file)
    merge!(S)
    detrend!(S)
    ungap!(S)
    taper!(S,t_max=50.)
    filtfilt!(S,rt="Bandpass",fl=freqmin,fh=freqmax,np=2)
    net,sta,loc,chan = split(S.id[1],'.')
    ind = findfirst(RESP.id .== S[1].id)
    S.loc[1] = RESP[ind].loc
    S.gain[1] = RESP[ind].gain
    S.resp[1] = RESP[ind].resp
    remove_resp!(S)
    return S
end

function all2all!(C::AbstractArray(CorrData),FFT1::AbstractArray{FFTData},FFT2::AbstractArray{FFTData},maxlag::Real)
    @assert size(FFT1,1) == 3
    @assert size(FFT2,1) == 3 
    @assert size(C,1) == 9 
    @inbounds for ii = 1:3
        @inbounds for jj = 1+3
            C[(ii-1) * 3 + jj] = correlate(FFT1[ii],FFT2[jj],maxlag)
        end
    end
    return nothing
end




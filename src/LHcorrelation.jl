export prepare_LH, LH_to_FFT, read_and_remove

function prepare_LH(fileE,fileN,fileZ,freqmin,freqmax,RESP)
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
    return E, N, Z 
end

function LH_to_FFT(E,N,Z,cc_len,cc_step,whitemin,whitemax)
    Eraw = RawData(E,cc_len,cc_step)
    Nraw = RawData(N,cc_len,cc_step)
    Zraw = RawData(Z,cc_len,cc_step)

    # simple processing 
    detrend!(Eraw)
    detrend!(Nraw)
    detrend!(Zraw)
    clip!(Eraw,3)
    clip!(Nraw,3)
    clip!(Zraw,3)

    # convert to fft 
    Efft = rfft(Eraw)
    Nfft = rfft(Nraw)
    Zfft = rfft(Zraw)

    # whiten E,N,Z components with 100 s amplitude smoothing 
    Efft.fft ./= smooth(abs.(Efft.fft),100)
    Nfft.fft ./= smooth(abs.(Efft.fft),100)
    Zfft.fft ./= smooth(abs.(Zfft.fft),100)
    return Efft, Nfft, Zfft 
end

function read_and_remove(file::String,freqmin::Real,freqmax::Real,RESP::SeisData)
    S = read_data("mseed",file)
    merge!(S)
    demean!(S)
    detrend!(S)
    ungap!(S)
    filtfilt!(S,rt="Bandpass",fl=freqmin,fh=freqmax,np=2)
    net,sta,loc,chan = split(S.id[1],'.')
    ind = findfirst(RESP.id .== S[1].id)
    S.loc[1] = RESP[ind].loc
    S.gain[1] = RESP[ind].gain
    S.resp[1] = RESP[ind].resp
    remove_resp!(S)
    return S
end



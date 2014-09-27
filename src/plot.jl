


# convert image to semilogy scaled image
function logimage(A,x)
    @assert size(A,1)==length(x)
    B = copy(A)
    y = logspace(log10(x[1]+0.1),log10(x[end]),length(x))
    for i=1:length(x)
        ind = searchsortedlast(x, y[i])
        ind < 1 && (ind=1)
        B[i,:] = A[ind,:]
    end
    return B
end

function logbinavg(p::AbstractArray, f::AbstractArray, bins::Integer)
    @assert length(p) == length(f)
    n = length(p)
    fls = logspace(log10(1), log10(f[end]), bins)
    x = zeros(bins)
    j = 1
    for i = 1:bins
        sum = 0
        c = 0
        while( j <= n && f[j] <= fls[i] )
            sum += p[j]
            c += 1
            j += 1
        end
        if c != 0
            x[i] = sum/c
        else
            x[i] = -200dB
        end
    end
    
    return x, fls
end

function plotall(s::MCSound; pgwin=tukey(samples(s),0.0003), pgnfft=nextfastfft(samples(s)), pgbins::Integer=600)

    fs = s.rate
    n = samples(s)
    x = s[:left]

    wl = 256*2
    wl > n && (wl=n>>2)
    gain = 2
    dbr = 70

    # audio
    sv = x
    st = [0:1/fs:(n-1)/fs]
    
    # periodogram
    P0 = 1 
    P = periodogram(x; onesided=true, nfft=pgnfft, fs=fs, window=pgwin)
    Pp, Pf = logbinavg(power(P), freq(P), pgbins)
    Pp = 10*log10(Pp./P0)
    
    #specrogram
    SG = spectrogram(x, wl, wl>>1; onesided=true, fs=fs, window=hanning)
    SGp = 10*log10(power(SG))
    SGp[find(x->(x<-dbr),SGp)] = -dbr
    SGp[find(x->(x>-gain),SGp)] = -gain
    SGt = time(SG)
    SGf = freq(SG)

    # PLOT

    f, ax = plt.subplots(3, 1, sharex=false)
    ax[1][:plot](st, sv, "b")
    ax[1][:set_xlim]([0,st[end]])
    ax[1][:set_ylim]([-1.1,1.1])
    #ax[1][:set_xlabel]("seconds")

    ax[2][:imshow](flipud(logimage(SGp,SGf)),aspect="auto",interpolation="bilinear",extent=[SGt[1],SGt[end],SGf[1]+0.1,SGf[end]])
    #ax[2][:pcolormesh](repmat(vec(SGt),1,length(SGf)),repmat(vec(SGf),length(SGt),1),flipud(SGp))
    ax[2][:set_ylim]([60,SGf[end]])
    ax[2][:set_yscale]("log")

    ax[2][:set_ylabel]("Hz")
    #ax[2][:set_xlabel]("seconds")

    ax[3][:semilogx](Pf, Pp, "b")
    ax[3][:set_xlim]([10,Pf[end]])
    ax[3][:set_ylim]([-100,2])
    #ax[3][:set_xlabel]("Hz")
    ax[3][:set_ylabel]("dB")

    #f[:set_size_inches](10,7)
    draw()

end

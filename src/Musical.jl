module Musical
export TIME_TYPE, AMP_TYPE, totimetype, toamptype,
    MCAudio, addto!, appendzeros!,
    nchannels, nsamples, 
    ismono, isstereo,
    left, right, 
    centerrange, leftrange, rightrange,
    seconds2saferange, unsafe_seconds2range,
    mono, stereo, 
    db, dB,
    importsound, exportsound,
    plotall,
    Envelope, ADSREnvelope, ENVELOPE_BOX, ENVELOPE_STANDARD, apply!, 
    sinewave, squarewave, sawtooth, trianglewave, whitenoise, binarynoise, chirpf, ichirpf,
    WaveForm, waveform,
    Note, Notes, scale!, render, 
    Tuning, TUNING_STANDARD, notefreq, bar
    
using WAV, DSP, PyPlot


typealias TIME_TYPE Float64
typealias AMP_TYPE Float32
totimetype(x) = convert(TIME_TYPE, x)
totimetype(x::Array) = convert(Array{TIME_TYPE}, x)
toamptype(x) = convert(AMP_TYPE, x)
toamptype(x::Array) = convert(Array{AMP_TYPE}, x)

# multi channel sound
type MCAudio
    v::Array{AMP_TYPE}
    fs::Int
    function MCAudio(v,fs)
        new(v,fs)
    end
end
function Base.getindex(s::MCAudio, op::Symbol)
    if op == :left
        @assert isstereo(s)
        return s.v[leftrange(s)]
    elseif op == :right
        @assert isstereo(s)
        return s.v[rightrange(s)]
    elseif op == :center
        @assert ismono(s)
        return s.v[centerrange(s)]
    else
        throw_unknown_symbol(op)
    end
end
function Base.getindex(s::MCAudio, op::Symbol, a::Real, b::Real)
    if op == :sec
        i1 = int(a*s.fs) + 1
        i1 = max(1, i1)
        i2 = int(b*s.fs)
        i2 = min(nsamples(s), i2)
        return s.v[i1:i2,:]
    else
        throw_unknown_symbol(op)
    end
end

function Base.setindex!(s::MCAudio, val::Real, op::Symbol)
    if op == :left
        @assert isstereo(s)
        range = leftrange(s)
    elseif op == :right
        @assert isstereo(s)
        range = rightrange(s)
    elseif op == :center
        @assert ismono(s)
        range = centerrange(s)
    else
        throw_unknown_symbol(op)
    end
    @inbounds for i in range
        s.v[i] = val
    end
    return s
end
function Base.setindex!(s::MCAudio, x::AbstractArray, op::Symbol)
    if nchannels(x) != 1 || nsamples(s) != nsamples(x)
        Base.throw_setindex_mismatch(x, (op,))
    end
    if op == :left
        @assert isstereo(s)
        si = leftrange(s)[1]
    elseif op == :right
        @assert isstereo(s)
        si = rightrange(s)[1]
    elseif op == :center
        @assert ismono(s)
        si = centerrange(s)[1]
    else
        throw_unknown_symbol(op)
    end
    copy!(s.v, si, x, 1, nsamples(x))
    return s
end

function Base.setindex!(s::MCAudio, val::Real, op::Symbol, a::Real, b::Real)
    if op == :sec
        @inbounds for j in 1:nchannels(s)
            for i in seconds2saferange(s, a, b)
                s.v[i,j] = val
            end
        end
        return s
    else
        throw_unknown_symbol(op)
    end
end
function Base.setindex!(s::MCAudio, x::AbstractArray, op::Symbol, a::Real, b::Real)
    if nchannels(s) != nchannels(x)
        Base.throw_setindex_mismatch(x, (op,))
    end
    if op == :sec
        r = seconds2saferange(s, a, b)
        if length(r) != nsamples(x)
            Base.throw_setindex_mismatch(x, (op,))
        end
        for j in 1:nchannels(s)
            copy!(s.v, r[1] + (j-1)*nsamples(s), x, 1 + (j-1)*nsamples(x), nsamples(x))
        end
        return s
    else
        throw_unknown_symbol(op)
    end
end

throw_unknown_symbol(x::Symbol) = error("uknown symbol :", string(x))

nchannels(s::MCAudio) = nchannels(s.v)
nchannels(v::AbstractArray) = size(v, 2)
nsamples(s::MCAudio) = nsamples(s.v)
nsamples(v::AbstractArray) = size(v, 1)
timelength(s::MCAudio) = totimetype(nsamples(s)/s.fs)
# mono
centerrange(s::MCAudio) = 1:nsamples(s)
ismono(s::MCAudio) = (nchannels(s) == 1)
# stereo
leftrange(s::MCAudio) = 1:nsamples(s)
rightrange(s::MCAudio) = nsamples(s)+1:2*nsamples(s)
isstereo(s::MCAudio) = (nchannels(s) == 2)

seconds2saferange(s::MCAudio, a::Real, b::Real) = seconds2saferange(s.fs,nsamples(s),a,b)
function seconds2saferange(fs::Integer, nsamples::Integer, a::Real, b::Real)
    i1 = int(totimetype(a)*fs) + 1
    i1 = max(1, i1)
    i2 = int(totimetype(b)*fs)
    i2 = min(nsamples, i2)
    return i1:i2
end
unsafe_seconds2range(s::MCAudio, a::Real, b::Real) = unsafe_seconds2range(s.fs,nsamples(s),a,b)
function unsafe_seconds2range(fs::Integer, nsamples::Integer, a::Real, b::Real)
    i1 = int(totimetype(a)*fs) + 1
    i2 = int(totimetype(b)*fs)
    return i1:i2
end

#Base.convert(::Type(MCAudio), x::Array) = MCAudio(x,fs)
function left(s::MCAudio) 
    @assert isstereo(s)
    return MCAudio(s[:left], s.fs)
end
function right(s::MCAudio) 
    @assert isstereo(s)
    return MCAudio(s[:right], s.fs)
end
function mono(s::MCAudio) 
    @assert nchannels(s) <= 2
    if ismono(s)
        return deepcopy(s)
    else
        return MCAudio((s[:left] + s[:right])/2, s.fs)
    end
end
function stereo(s::MCAudio) 
    @assert nchannels(s) <= 2
    if ismono(s)
        return MCAudio(repmat(s.v,1,2), s.fs)
    else
        return deepcopy(s)
    end
end
function stereo(s1::MCAudio, s2::MCAudio) 
    @assert ismono(s1) && ismono(s2)
    @assert s1.fs == s2.fs
    return MCAudio([s1.v s2.v], s1.fs)
end

# + - * / ^
+(s::MCAudio, a::Union(Real, AbstractArray)) = MCAudio(s.v .+ a, s.fs)
+(a::Union(Real, AbstractArray), s::MCAudio) = MCAudio(a .+ s.v, s.fs)
-(s::MCAudio, a::Union(Real, AbstractArray)) = MCAudio(s.v .- a, s.fs)
-(a::Union(Real, AbstractArray), s::MCAudio) = MCAudio(a .- s.v, s.fs)
*(s::MCAudio, a::Union(Real, AbstractArray)) = MCAudio(s.v .* a, s.fs)
*(a::Union(Real, AbstractArray), s::MCAudio) = MCAudio(a .* s.v, s.fs)
/(s::MCAudio, a::Union(Real, AbstractArray)) = MCAudio(s.v ./ a, s.fs)
/(a::Union(Real, AbstractArray), s::MCAudio) = MCAudio(a ./ s.v, s.fs)
#(^)(s::MCAudio, a::Real) = MCAudio(s.v .^ a, s.fs) #defined in Base
(^)(a::Union(Real, AbstractArray), s::MCAudio) = MCAudio(a .^ s.v, s.fs)
+(s::MCAudio, d::MCAudio) = MCAudio(s.v .+ d.v, s.fs)
-(s::MCAudio, d::MCAudio) = MCAudio(s.v .- d.v, s.fs)
*(s::MCAudio, d::MCAudio) = MCAudio(s.v .* d.v, s.fs)
/(s::MCAudio, d::MCAudio) = MCAudio(s.v ./ d.v, s.fs)
(^)(s::MCAudio, d::MCAudio) = MCAudio(s.v .^ d.v, s.fs)


function addto!(s::MCAudio, time::Real, a::MCAudio)
    @assert nchannels(s) == nchannels(a)
    @assert time >= 0

    starti = unsafe_seconds2range(s, time, time).start
    ir = starti:(starti + nsamples(a) - 1)
    if ir[end] > nsamples(s)  #append s with zeros
        appendzeros!(s, ir[end])
    end
    addarrays!(s.v, ir, a.v)
    return s
end
function appendzeros!(s::MCAudio, n::Integer)
    @assert n > nsamples(s)
    sv = s.v
    s.v = zeros(AMP_TYPE, n, nchannels(s))
    for j in 1:nchannels(s)
        copy!(s.v, 1 + (j-1)*nsamples(s.v), sv, 1 + (j-1)*nsamples(sv), nsamples(sv))
    end
    return s
end
function addarrays!(a::Array, ir::Range, b::Array)
    @assert nchannels(a) == nchannels(b)
    @assert length(ir) == nsamples(b)
    @assert ir[end] <= nsamples(a)
    @inbounds for c in 1:nchannels(a)
        for i = 1:length(ir)
            a[ir[i],c] += b[i,c]
        end
    end
    return a
end


# utils

type dBconvert
end
const dB = dBconvert()
*(a::Real, d::dBconvert) = db(a)
db(a::Real) = toamptype(10^(a*0.05))


function importsound(file::String; subrange=Any)
    x = wavread(file; subrange=subrange, format="double")
    return MCAudio(toamptype(x[1]), convert(Int,x[2]))
end
function exportsound(s::MCAudio, file::String)
    wavwrite(s.v, file; Fs=s.fs, nbits=32, compression=WAVE_FORMAT_IEEE_FLOAT)
end

include("plot.jl")

include("synth.jl")



end #module


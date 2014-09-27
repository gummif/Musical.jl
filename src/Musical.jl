module Musical
export MCSound, 
    nchannels, samples, 
    ismono, isstereo,
    left, right, 
    centerrange, leftrange, rightrange,
    seconds2saferange,
    mono, stereo, 
    db, dB,
    importsound, exportsound,
    Envelope, ADSREnvelope, ENVELOPE_BOX, ENVELOPE_STANDARD, apply!, 
    sinewave, squarewave, sawtooth, trianglewave, whitenoise, binarynoise, waveform,
    plotall
using WAV, DSP, PyPlot
require("Options")
using OptionsMod


# multi channel sound
type MCSound
    v::Array{Float32}
    rate::Integer
    function MCSound(v,rate)
        new(v,rate)
    end
end
function Base.getindex(s::MCSound, op::Symbol)
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
        error("unknown symbol")
    end
end
function Base.getindex(s::MCSound, op::Symbol, a::Real, b::Real)
    if op == :sec
        i1 = int(a*s.rate) + 1
        i1 = max(1, i1)
        i2 = int(b*s.rate)
        i2 = min(samples(s), i2)
        return s.v[i1:i2,:]
    else
        error("unknown symbol")
    end
end

function Base.setindex!(s::MCSound, val::Real, op::Symbol)
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
        error("unknown symbol")
    end
    @inbounds for i in range
        s.v[i] = val
    end
    return s
end
function Base.setindex!(s::MCSound, x::AbstractArray, op::Symbol)
    if nchannels(x) != 1 || samples(s) != samples(x)
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
        error("unknown symbol")
    end
    copy!(s.v, si, x, 1, samples(x))
    return s
end

function Base.setindex!(s::MCSound, val::Real, op::Symbol, a::Real, b::Real)
    if op == :sec
        @inbounds for j in 1:nchannels(s)
            for i in seconds2saferange(s, a, b)
                s.v[i,j] = val
            end
        end
        return s
    else
        error("unknown symbol")
    end
end
function Base.setindex!(s::MCSound, x::AbstractArray, op::Symbol, a::Real, b::Real)
    if nchannels(s) != nchannels(x)
        Base.throw_setindex_mismatch(x, (op,))
    end
    if op == :sec
        r = seconds2saferange(s, a, b)
        if length(r) != samples(x)
            Base.throw_setindex_mismatch(x, (op,))
        end
        for j in 1:nchannels(s)
            copy!(s.v, r[1] + (j-1)*samples(s), x, 1 + (j-1)*samples(x), samples(x))
        end
        return s
    else
        error("unknown symbol")
    end
end

nchannels(s::MCSound) = nchannels(s.v)
nchannels(v::AbstractArray) = size(v, 2)
samples(s::MCSound) = samples(s.v)
samples(v::AbstractArray) = size(v, 1)
# mono
centerrange(s::MCSound) = 1:samples(s)
ismono(s::MCSound) = (nchannels(s) == 1)
# stereo
leftrange(s::MCSound) = 1:samples(s)
rightrange(s::MCSound) = samples(s)+1:2*samples(s)
isstereo(s::MCSound) = (nchannels(s) == 2)

seconds2saferange(s::MCSound, a::Real, b::Real) = seconds2saferange(s.rate,samples(s),a,b)
function seconds2saferange(rate::Integer, samples::Integer, a::Real, b::Real)
    i1 = int(a*rate) + 1
    i1 = max(1, i1)
    i2 = int(b*rate)
    i2 = min(samples, i2)
    return i1:i2
end

#Base.convert(::Type(MCSound), x::Array) = MCSound(x,rate)
function left(s::MCSound) 
    @assert isstereo(s)
    return MCSound(s[:left], s.rate)
end
function right(s::MCSound) 
    @assert isstereo(s)
    return MCSound(s[:right], s.rate)
end
function mono(s::MCSound) 
    @assert nchannels(s) <= 2
    if ismono(s)
        return deepcopy(s)
    else
        return MCSound((s[:left] + s[:right])/2, s.rate)
    end
end
function stereo(s::MCSound) 
    @assert nchannels(s) <= 2
    if ismono(s)
        return MCSound(repmat(s.v,1,2), s.rate)
    else
        return deepcopy(s)
    end
end
function stereo(s1::MCSound, s2::MCSound) 
    @assert ismono(s1) && ismono(s2)
    @assert s1.rate == s2.rate
    return MCSound([s1.v s2.v], s1.rate)
end

# + - * / ^
+(s::MCSound, a::Real) = MCSound(s.v .+ a, s.rate)
+(a::Real, s::MCSound) = MCSound(a .+ s.v, s.rate)
-(s::MCSound, a::Real) = MCSound(s.v .- a, s.rate)
-(a::Real, s::MCSound) = MCSound(a .- s.v, s.rate)
*(s::MCSound, a::Real) = MCSound(s.v .* a, s.rate)
*(a::Real, s::MCSound) = MCSound(a .* s.v, s.rate)
/(s::MCSound, a::Real) = MCSound(s.v ./ a, s.rate)
/(a::Real, s::MCSound) = MCSound(a ./ s.v, s.rate)
#(^)(s::MCSound, a::Real) = MCSound(s.v .^ a, s.rate) #defined in Base
(^)(a::Real, s::MCSound) = MCSound(a .^ s.v, s.rate)
+(s::MCSound, d::MCSound) = MCSound(s.v .+ d.v, s.rate)
-(s::MCSound, d::MCSound) = MCSound(s.v .- d.v, s.rate)
*(s::MCSound, d::MCSound) = MCSound(s.v .* d.v, s.rate)
/(s::MCSound, d::MCSound) = MCSound(s.v ./ d.v, s.rate)
(^)(s::MCSound, d::MCSound) = MCSound(s.v .^ d.v, s.rate)

type dBconvert
end
const dB = dBconvert()
*(a::Real, d::dBconvert) = db(a)
db(a::Real) = float32(10^(a*0.05))

function importsound(file::String; subrange=Any)
    x = wavread(file; subrange=subrange, format="double")
    return MCSound(float32(x[1]), convert(Int,x[2]))
end
function exportsound(s::MCSound, file::String)
    wavwrite(s.v, file; Fs=s.rate, nbits=32, compression=WAVE_FORMAT_IEEE_FLOAT)
end

include("plot.jl")



# Synthesizer

abstract Envelope
type ADSREnvelope <: Envelope
    A::Float32      #time
    D::Float32      #time
    S::Float32      #amp
    R::Float32      #time
    amp::Float32    #attack
    ADSREnvelope(A,D,S,R,amp) = new(A,D,S,R,amp)
end
const ENVELOPE_BOX = ADSREnvelope(0,0,1,0,1)
const ENVELOPE_STANDARD = ADSREnvelope(3f-2,1f-1,7f-1,5f-1,9f-1)

releasetime(env::ADSREnvelope) = env.R
function apply!(w::Array, env::ADSREnvelope, hold::Real, fs::Integer)
    hold = float32(hold)
    t1 = env.A
    t2 = env.A + env.D
    t3 = hold
    t4 = hold + env.R
    if t1 < hold < t2
        t2 = hold
    elseif hold < env.A
        t1 = hold
        t2 = hold
    end
    r = seconds2saferange(fs,samples(w),0f0,t1)
    w[r,:] .*= linspace(0f0,env.amp,length(r))
    r = seconds2saferange(fs,samples(w),t1,t2)
    w[r,:] .*= linspace(env.amp,env.S,length(r))
    r = seconds2saferange(fs,samples(w),t2,t3)
    w[r,:] .*= linspace(env.S,env.S,length(r))
    r = seconds2saferange(fs,samples(w),t3,t4)
    w[r,:] .*= linspace(env.S,0f0,length(r))

    return w
end



# waves are periodic functions with period 1, and : [0, inf) -> [-1,1]
const twopi = float32(2*pi)
const twodpi = float32(2/pi)
nozerosign(x) = ifelse(x >= 0, one(x), oftype(x,-1))

sinewave(x) = sin(twopi*x)
squarewave(x) = nozerosign(sin(twopi*x))
sawtooth(x) = 2*(x + 5f-1 - floor(x + 5f-1)) - 1
trianglewave(x) = twodpi*asin(sin(twopi*x))
whitenoise(x) = 2*rand(Float32) - 1
binarynoise(x) = 2*round(rand(Float32)) - 1


function float32wave(wave::Function, r::Range, channels::String="Mono", imaging::Real=0)
    if channels == "Mono"
        w = Array(Float32, length(r))
        c = 1
        for i in r
            w[c] = wave(i)
            c += 1
        end
    elseif channels == "Stereo"
        @assert imaging >= 0
        w = Array(Float32, length(r), 2)
        c = 1
        for i in r
            w[c,1] = wave(i)
            c += 1
        end
        c = 1
        for i in r+imaging
            w[c,2] = wave(i)
            c += 1
        end
    else
        error("unknown channels")
    end
    return w
end

# if wave is a periodic function with period 1
# returns enveloped waveform with frequency freq
function waveform(freq::Real, time::Real, fs::Integer, opts::Options=@options __=nothing)
    @defaults opts wave=sinewave env=ENVELOPE_BOX channels="Mono" imaging=0
    wave::Function
    env::Envelope
    channels::String
    imaging::Real
    r = float32((0:(fs-1)*(time+releasetime(env)))*freq/fs)
    w = float32wave(wave, r, channels, imaging)
    apply!(w, env, time, fs)
    return MCSound(w, fs)
end



end


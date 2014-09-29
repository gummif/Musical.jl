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
    Note, Notes, Tuning, TUNING_STANDARD, notefreq
    
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



# Synthesizer

abstract Envelope
type ADSREnvelope <: Envelope
    A::TIME_TYPE      #time
    D::TIME_TYPE      #time
    S::AMP_TYPE       #amp
    R::TIME_TYPE      #time
    amp::AMP_TYPE     #attack
    ADSREnvelope(A,D,S,R,amp) = new(A,D,S,R,amp)
end
const ENVELOPE_BOX = ADSREnvelope(0,0,1,0,1)
const ENVELOPE_STANDARD = ADSREnvelope(3f-2,1f-1,7f-1,5f-1,9f-1)
releasetime(env::ADSREnvelope) = env.R

function apply!(w::Array, env::ADSREnvelope, hold::Real, fs::Integer)
    hold = totimetype(hold)
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
    applylinspace!(w, fs, nsamples(w),  0, t1, toamptype(0), env.amp)
    applylinspace!(w, fs, nsamples(w), t1, t2, env.amp, env.S)
    applylinspace!(w, fs, nsamples(w), t2, t3, env.S, env.S)
    applylinspace!(w, fs, nsamples(w), t3, t4, env.S, toamptype(0))
    return w
end
function applylinspace!(w, fs, n, t1, t2, A1, A2)
    r = seconds2saferange(fs, n, t1, t2)
    m = ifelse( length(r)==1, A2, (A2 - A1)/(length(r) - 1) )
    @inbounds for c in 1:nchannels(w)
        for i = 1:length(r)
            w[r[i],c] *= A1 + (i-1)*m
        end
    end
end


# waves are periodic functions with period 1, and : [0, inf) -> [-1,1]
const twopi = totimetype(2*pi)
const twodpi = totimetype(2/pi)
const halft = totimetype(0.5)
nozerosign(x) = ifelse(x >= 0, one(x), oftype(x,-1))

sinewave(x) = sin(twopi*x)
squarewave(x) = nozerosign(sin(twopi*x))
sawtooth(x) = 2*(x + halft - floor(x + halft)) - 1
trianglewave(x) = twodpi*asin(sin(twopi*x))
whitenoise(x) = 2*rand(AMP_TYPE) - 1
binarynoise(x) = 2*round(rand(AMP_TYPE)) - 1
chirpf(x, k::Real=1, f0::Real = 0) = x.*(f0 + k*x)
ichirpf(x) = chirpf(x, -1, 20000)
expchirpf(x, k::Real=2, f0::Real = 1) = x.*(1f-2.*k.^(x))


type WaveForm
    wave::Function              # wave function
    env::Envelope               # envelope
    nchannels::Int              # number of channels
    phase::Vector{TIME_TYPE}    # phase shift of each channel
    function WaveForm(wave ,env, nchannels, phase)
        isa(phase, Real) && ( phase = repmat([phase], nchannels) )
        @assert nchannels == length(phase)
        return new(wave, env, nchannels, phase)
    end
end
WaveForm() = WaveForm(sinewave, ENVELOPE_BOX, 1, [0])
function WaveForm(; wave::Function=sinewave, 
                env::Envelope=ENVELOPE_BOX, 
                nchannels::Integer=1, 
                phase::Union(Real,Vector)=0)
    WaveForm(wave, env, nchannels, phase)
end

# if wave is a periodic function with period 1
# returns enveloped waveform with frequency freq
function waveform(freq::Real, time::Real, fs::Integer, form::WaveForm=WaveForm())
    r = (0:(fs-1)*(totimetype(time)+releasetime(form.env)))*totimetype(freq/fs)
    w = time2wave(r, form)
    apply!(w, form.env, time, fs)
    return MCAudio(w, fs)
end

# return an AMP_TYPE Array
function time2wave(r::Range, form::WaveForm)
    n = length(r)
    if form.nchannels == 1
        w = Array(AMP_TYPE, n)
    else
        w = Array(AMP_TYPE, n, form.nchannels)
    end
    ww = Array(AMP_TYPE, n)
    return setwave!(w, ww, r, form.wave, form.phase)
end
function setwave!(w, ww, r, wave, fp)
    n = length(r)
    for c in 1:nchannels(w)
        rp = r + fp[c]
        broadcast!(wave, ww, rp)        # ww = wave(rp)
        copy!(w, 1 + (c-1)*n, ww, 1, n) # w[:,c] = ww
    end
    return w
end

type Note
	freq::Real
	start::TIME_TYPE	# starting time
	S::TIME_TYPE		# sustain time
	amp::AMP_TYPE		# amplitude multiplier
	Note(freq, start, S, amp=1) = new(freq, start, S, amp)
end

type Notes
	notes::Vector{Note}
	wf::WaveForm
	fs::Int
    Notes(notes, wf, fs) = new(notes, wf, fs)
end

immutable Tuning
    value2freq::Function    # convert frequency value to Hz (Real -> TIME_TYPE)
    note2value::Function    # convert note name to frequency value (ASCIIString -> Real)
    Tuning(value2freq, note2value) = new(value2freq, note2value)
end

function midinotenumber2freq(d::Real)
    return totimetype( 2^((d-69)/12)*440 )
end
function standardnote2value(n::ASCIIString)
    v = 60
    baseoct = 4
    note, oct = noteoct(n, baseoct)
    v = v + (oct - baseoct)*12
    if note == "C"
        # v += 0
    elseif note == "C#" || note == "Db"
        v += 1 
    elseif note == "D"
        v += 2
    elseif note == "D#" || note == "Eb"
        v += 3
    elseif note == "E"
        v += 4
    elseif note == "F"
        v += 5
    elseif note == "F#" || note == "Gb"
        v += 6
    elseif note == "G"
        v += 7
    elseif note == "G#" || note == "Ab"
        v += 8
    elseif note == "A"
        v += 9
    elseif note == "A#" || note == "Bb"
        v += 10
    elseif note == "B"
        v += 11
    else
        error("uknown note ", string(note))
    end
    return v
end
function noteoct(n::ASCIIString, baseoct::Integer)
    ind = 0
    for i in 1:length(n)
        if isdigit(n[i])
            ind = i
            break
        end
    end
    if ind == 0
        oct = baseoct
    elseif ind == 1
        error("no note supplied: ", string(n))
    else
        oct = int(n[ind:end])
    end
    note = n[1:ind-1]
    return note, oct
end

const TUNING_STANDARD = Tuning(midinotenumber2freq, standardnote2value)

notefreq(x::ASCIIString, t::Tuning=TUNING_STANDARD) = notefreq(t.note2value(x), t)
notefreq(x::Real, t::Tuning=TUNING_STANDARD) = t.value2freq(x)



end #module


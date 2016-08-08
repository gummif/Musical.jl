

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
    @assert time >= 0
    r = (0:(fs-1)*(totimetype(time)+releasetime(form.env)))*totimetype(freq/fs)
    w = time2wave(r, form)
    apply!(w, form.env, time, fs)
    return MCAudio(w, fs)
end
waveform(tone::ASCIIString, args...) = waveform(notefreq(tone), args...)


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
	freq::TIME_TYPE     # frequency
	start::TIME_TYPE	# starting time
	S::TIME_TYPE		# sustain time
	amp::AMP_TYPE		# amplitude multiplier
	Note(freq, start, S, amp=1) = new(freq, start, S, amp)
end
Note(tone::ASCIIString, args...) = Note(notefreq(tone), args...)

type Notes
	notes::Vector{Note}
	wf::WaveForm
	fs::Int
    Notes(notes, wf, fs) = new(notes, wf, fs)
end

function Base.scale!(n::Note, from::Real, to::Real)
    n.freq *= to/from 
    return n
end
function Base.scale!(ns::Notes, from::Real, to::Real)
    for i in 1:length(ns.notes)
        scale!(ns.notes[i], from, to)
    end
    return ns
end

# convert Notes to MCAudio
function render(ns::Notes)
    maxt = findmaxtime(ns.notes)
    r = unsafe_seconds2range(ns.fs, 0, 0.0, maxt)
    fs = ns.fs
    wf = ns.wf
    s = MCAudio(zeros(AMP_TYPE, r[end], ns.wf.nchannels), fs)

    for i in 1:length(ns.notes)
        w = waveform(ns.notes[i].freq, ns.notes[i].S, fs, wf)
        scale!(ns.notes[i].amp, w.v)
        addto!(s, ns.notes[i].start, w)
    end
    return s
end
function findmaxtime(ns::Vector{Note})
    maxt = 0
    @inbounds for i in 1:length(ns)
        @assert ns[i].start >= 0
        ti = ns[i].start + ns[i].S
        ti > maxt && (maxt = ti)
    end
    return maxt
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
        error("unknown note ", string(note))
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

# convert bar to seconds
function bar(bt::Tuple, bpm::Real=120, sig::NTuple{2,Integer}=(4,4))
    n = length(bt)
    if n > 1
        for i in 2:n
            @assert 1 <= bt[i] < sig[1] + 1
        end
    end
    @assert sig == (4,4)  # other not implemented
    spb = 60/totimetype(bpm)    # seconds per beat
    t = totimetype(0)
    if n >= 1
        t += (totimetype(bt[1])-1)*totimetype(sig[1])*spb
    end
    if n >= 2
        t += (totimetype(bt[2])-1)*spb
    end
    if n >= 3
        t += (totimetype(bt[3])-1)*spb/totimetype(sig[1])
    end
    n >= 4 && error("tuple to long")
    
    return t::TIME_TYPE
end





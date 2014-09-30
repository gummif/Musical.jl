# Musical.jl

Manipulate multi-channel audio in an intuitive way using the programming power of Julia. Requires the Julia packages WAV, DSP, and PyPlot.

All audio is stored as Float32 and the default export format is the IEEE 32bit float. For time, Float64 is used internally since Float32 is very inaccurate above a few hundred seconds.

For sample-rate and bit-rate convertions it is advised to use other software.

See license (MIT) in LICENSE.md.


## Usage

```julia
# import audio as a MCAudio object containing Float32 data
s = importsound("sounds.wav")
s.v     # the array of values
s.fs    # the sample-rate
# extract the left channel of stereo audio
lc = s[:left]
# divide the first 100ms by 2
s[:sec, 0, 0.1] /= 2
# multiply the right channel of s by the squared left channel
s[:right] .*= s[:left].^2
# create a 0.3s sustained 200Hz waveform at 44100Hz sample-rate,
# with a standard ADSR envelope, decreased by -3dB, with 0.2*wavelength stereo imaging
wf = WaveForm(wave=sawtooth, env=ENVELOPE_STANDARD, nchannels=2, phase=[0, 0.2])
w = -3dB*waveform(200, 0.3, 44100, wf)
# plotting utilities
plotall(w)
# export
exportsound(w,"waveform.wav")
```

Example of creation of notes in standard tuning
```julia
wf = WaveForm(env=ENVELOPE_STANDARD)
fs = 44100
w = waveform(notefreq("E4"), 5, fs, wf)+0.5*waveform(notefreq("E5"), 5, fs, wf)+0.3*waveform(notefreq("E6"), 5, fs, wf)
w = w/2
exportsound(w, "E4.wav")
```

Example of creation of a sequence of notes in standard tuning and plotting
```julia
using Musical
env = ADSREnvelope(0.07,0.05,0.8,0.8,1)
wf = WaveForm(wave=x->trianglewave(x), env=env, nchannels=2)
fs = 44100
S = 0.7
bpm = 150
notevec = [ Note("A3", bar((1,1),bpm), S), 
			Note("A4", bar((2,1),bpm), S), 
			Note("C4", bar((3,1),bpm), S), 
			Note("F#4", bar((4,1),bpm), S)]
ns = Notes(notevec, wf ,fs)
s = render(ns)
ss = render(scale!(ns, notefreq("G3"), notefreq("F3")))
addto!(s, bar((5,), bpm), ss)

exportsound(s,"notes.wav")

f,ax = plotall(s)

using PyPlot
f[:set_size_inches](6,5)
savefig("plot.png")
```

![Plot](/test/plot.png)





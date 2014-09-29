

using Musical, DSP
sawwave(x) = x - floor(x)
mywave(x) = 0.5*sin(2*pi*x) + 0.5*sign(sin(2*pi*x*2)) + 0.2*sin(2*pi*x*4+1)+ 0.2*sin(2*pi*x*4.02+1)
#mywave(x) = 1.1*sin(2*pi*x) + 0.3*sign(sin(2*pi*x))+sawwave(x*2)*0.7
expchirpf(x, k::Real=1.1, f0::Real = 1f-1) = x.*(f0.*k.^(x))
mywave(x) = sinewave(expchirpf(x))
env = ADSREnvelope(0.05,0.05,0.9,1,1)
wf = WaveForm(wave=x->sinewave(expchirpf(x)), env=env, nchannels=2, phase=[0, 0.2])

w = waveform(16, 7, 44100, wf)

w2 = waveform(16, 1, 44100, wf)
addto!(w, 6.5, w2)

exportsound(w,"test44b.wav")
#w = importsound("asmz.wav")


run(`vlc test44b.wav`)
#plotall(stereo(w))
plotall(addto!(w, 0, (-6dB)*waveform(2000, 1, 44100, WaveForm(nchannels=2))))

#f, ax = plotall(stereo(w))
#ax[3][:grid](axis="y", which="major")



wf = WaveForm(env=ENVELOPE_STANDARD)
w = waveform(notefreq("E4"), 5, 44100, wf)+0.5*waveform(notefreq("E5"), 5, 44100, wf)+0.3*waveform(notefreq("E6"), 5, 44100, wf)
w = w/2
w = waveform(notefreq("E2"), 5, 44100, wf)+0.5*waveform(notefreq("E3"), 5, 44100, wf)+0.3*waveform(notefreq("E4"), 5, 44100, wf)
w = w/2
w = waveform(notefreq("A2"), 5, 44100, wf)+0.5*waveform(notefreq("A3"), 5, 44100, wf)+0.3*waveform(notefreq("A4"), 5, 44100, wf)
w = w/2
w = waveform(notefreq("D3"), 5, 44100, wf)+0.5*waveform(notefreq("D4"), 5, 44100, wf)+0.3*waveform(notefreq("D5"), 5, 44100, wf)
w = w/2
w = waveform(notefreq("G3"), 5, 44100, wf)+0.5*waveform(notefreq("G4"), 5, 44100, wf)+0.3*waveform(notefreq("G5"), 5, 44100, wf)
w = w/2
w = waveform(notefreq("B3"), 5, 44100, wf)+0.5*waveform(notefreq("B4"), 5, 44100, wf)+0.3*waveform(notefreq("B5"), 5, 44100, wf)
w = w/2
exportsound(w, "bobo.wav")
wf = WaveForm(wave=sawtooth, env=ENVELOPE_STANDARD, nchannels=2, phase=[0, 0.1])
w = waveform(notefreq("G2"), 5, 44100, wf)+1*waveform(1.01*notefreq("D3"), 5, 44100, wf)+1.3*waveform(notefreq("G3"), 5, 44100, wf)+1*waveform(1.05*notefreq("A#3"), 5, 44100, wf)
w = w
exportsound(w, "bobo.wav")








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

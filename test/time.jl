
using Musical, PyPlot
env = ADSREnvelope(0.05,0.05,0.9,1,1)
wf = WaveForm(wave=sinewave, env=env, nchannels=2, phase=[0, 0.2])
w = waveform(200, 1, 44100, wf)

f() = for i = 1:100; waveform(200, 1, 44100, wf); end
f();
@time f();

@profile f();
Profile.print()

plot(w.v)



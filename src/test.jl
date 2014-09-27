

using Musical, DSP
sawwave(x) = x - floor(x)
mywave(x) = 0.5*sin(2*pi*x) + 0.5*sign(sin(2*pi*x*2)) + 0.2*sin(2*pi*x*4+1)+ 0.2*sin(2*pi*x*4.02+1)
#mywave(x) = 1.1*sin(2*pi*x) + 0.3*sign(sin(2*pi*x))+sawwave(x*2)*0.7
mywave(x) = sinewave(x)
env = ADSREnvelope(0.05,0.05,0.9,1,1)
wf = WaveForm(wave=mywave, env=env, nchannels=2, phase=[0, 0.2])
w = waveform(200, 1, 44100, wf)

exportsound(w,"test44b.wav")
#w = importsound("asmz.wav")

#using PyPlot
#plot(w[:left])

Fs = 44100;   
t = 0:1/Fs:.3;
x = cos(2*pi*t*14000)+0.1*randn(length(t)); 
#w = MCAudio(x,Fs)

plotall(stereo(w))



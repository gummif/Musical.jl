# Musical.jl

Manipulate multi-channel audio in an intuitive way using the programming power of Julia. Requires the Julia packages WAV, DSP, and PyPlot.

All audio is stored as Float32 and the default export format is the IEEE 32bit float. For samplingrate and bitrate changes it is advised to use other software.



## Usage

```julia
# import audio as a MCSound object containing Float32 data
s = importsound("sounds.wav")
s.v     # the array of values
s.rate  # the samplerate
# extact left channels of stereo audio
lc = s[:left]
# divide the first 100ms by 2
s[:sec, 0, 0.1] /= 2
# multiply the right channel of s by the squared left channel
s[:right] .*= s[:left].^2
# create a 0.3s sustained 200Hz waveform at 44100Hz samplerate,
# with a standard ADSR envelope, decreased by -3dB, with 0.2*wavelength stereo imaging
using OptionsMod
opts = @options wave=sawtooth env=ENVELOPE_STANDARD channels="Stereo" imaging=0.2
w = -3dB*waveform(200, 0.3, 44100, opts)
# plotting utilities
plotall(w)
# export
exportsound(w,"sounds2.wav")
```



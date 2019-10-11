
function read_as_axisarray(stream::LibSndFile.SndFileSource, seconds)
    s = SampledSignals.s
    Fs = SampledSignals.samplerate(stream)
    x = read(stream, seconds * s)
    tt = SampledSignals.domain(x)
    x = convert.(Float32, vec(x.data))
    x = AxisArray(x, Axis{:time}(tt))
    return x
end

function read_as_axisarray(stream::LibSndFile.SndFileSource)
    return read_as_axisarray(stream, seconds)
end

function read_as_axisarray(wavfile::AbstractString)
    stream = loadstreaming(wavfile)
    seconds = SampledSignals.nframes(stream) / SampledSignals.samplerate(stream)
    return read_as_axisarray(stream, seconds)
end

function read_as_axisarray(wavfile::AbstractString, seconds)
    stream = loadstreaming(wavfile)
    return read_as_axisarray(stream, seconds)
end

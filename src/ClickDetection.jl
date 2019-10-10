module ClickDetection


using AxisArrays
using IntervalSets
using DSP
# using MelGeneralizedCepstrums
using FileIO
import LibSndFile
import SampledSignals
import DataFrames
using Dates
using Statistics
using Plots

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
    seconds = SampledSignals.nframes(stream) / SampledSignals.samplerate(stream)
    return read_as_axisarray(stream, seconds)
end


struct ClickPointer
    series::AxisArray
    interval::ClosedInterval
    threshold::Real
end

left(c::ClickPointer) = c.interval.left
right(c::ClickPointer) = c.interval.right

function moving_avg(x, k::Integer)
    @assert isodd(k)
    n = length(x)
    x_smooth = zeros(n)
    w = floor(Int, k/2)
    for i in 1:n
        lower = max(1, i - w)
        upper = min(n, i + w)
        x_smooth[i] = mean(x[lower:upper])
    end
    return x_smooth
end

function teager(x::AbstractVector)
    @assert length(x) >= 3 "Not enough measurements to compute Teager energy"
    output = zeros(size(x))
    for i in 2:length(x)-1
        output[i] = x[i]^2 - x[i-1] * x[i+1]
    end
    output[1] = output[2]
    output[end] = output[end-1]
    return output
end

"""
Combine intervals which overlap.  Takes a vector of ClosedIntervals
"""
function merge_intervals(ints::Vector{I}, sort_first=true) where {I <: ClosedInterval}
    if length(ints) < 1
        return I[]
    end

    if sort_first
        sorted_ints = sort(ints, by = x -> x.left)
    end
    current_interval = sorted_ints[1]
    merged_intervals = I[]
    i = 2
    while i <= length(sorted_ints)
        if isempty(intersect(current_interval, sorted_ints[i]))
            push!(merged_intervals, current_interval)
            current_interval = sorted_ints[i]
        else
            current_interval = union(current_interval, sorted_ints[i])
        end
        i += 1
    end
    push!(merged_intervals, current_interval)
    return merged_intervals
end

function detect_clicks(x, thresh, span)
    x_teager = teager(x)
    times = axisvalues(x[Axis{:time}])[1]
    t0, t1 = times[1], times[end]
    high_teager_idx = findall(x -> x > thresh, x_teager)
    intervals = [max(t0, times[i]-span)..min(times[i]+span, t1) for i in high_teager_idx]
    intervals = merge_intervals(intervals)
    return [ClickPointer(x, i, thresh) for i in intervals]
end


function get_data(click::ClickPointer)
    x = click.series[click.interval]
    t = axisvalues(x)[1]
    return t, x
end

function get_background(click::ClickPointer, width)
    left = click.interval.left
    right = click.interval.right
    x = click.series
    background = [x[(left-width)..left]; x[right..(right+width)]]
    return background
end


function peak_lead_ratio(click::ClickPointer, nsamples_lead)
    t, x = get_data(click)
    x² = abs2.(x)
    return mean(x²[nsamples_lead:end]) / mean(x²[1:nsamples_lead])
end

function summarize_click(click::ClickDetector.ClickPointer, start_time, Fs)
     t, x = get_data(click)
     n = div(length(x), 8)
     pg = welch_pgram(x, n, fs=Fs)
     peak_pres = maximum(abs.(x))
     peak_freq = pg.freq[findmax(pg.power)[2]]
     mean_freq = sum(pg.freq .* pg.power / sum(pg.power))
     duration = click.interval.right - click.interval.left
     nzeros = zero_crossings(x)
     cepstrum = Array(estimate(LinearCepstrum(14), x))'
     datetime = start_time + Dates.Millisecond(round(click.interval.left * 1e3))
     return [datetime duration peak_pres peak_freq mean_freq nzeros cepstrum...]
end


function process_wavfile(wf; noise_pct=0.6, noise_factor=50, window=1e-3, lead_snr=5)
     start_time = DateTime(basename(wf)[6:20], "yyyymmdd_HHMMSS")
     println("Reading $wf...")
     wav_stream = loadstreaming(wf)
     Fs = SampledSignals.samplerate(wav_stream)
     x = read_as_axisarray(wav_stream)
     tt = axisvalues(x)[1]

     # High pass filter
     butter = digitalfilter(Highpass(5e3, fs=Fs), Butterworth(3))
     x = filtfilt(butter, x)

     # Notch filter to remove 50 kHz hum
     notch = iirnotch(50e3, 1e3, fs=Fs)
     x = filt(notch, x)
     x = AxisArray(x, Axis{:time}(tt))

     # Teager-engerygy click detection
     x_teager = teager(x)
     noise_floor = quantile(x_teager, noise_pct)
     thresh = noise_floor * noise_factor
     clicks = detect_clicks(x, thresh, window)

     # eliminate "clicks" where high-teager detection was just a random spike
     nsamples_lead = floor(Int, window * Fs)
     clicks = filter(c -> peak_lead_ratio(c, nsamples_lead) > lead_snr, clicks)

     # DataFrame for output
     starts = [start_time + Dates.Millisecond(round(left(c) * 1e3)) for c in clicks]
     stops = [start_time + Dates.Millisecond(round(right(c) * 1e3)) for c in clicks]
     waveforms = [get_data(c)[2] for c in clicks]
     if length(clicks) == 0
         return DataFrames.DataFrame(wavfile = wf, start = DataFrames.missing,
            stop = DataFrames.missing, waveform = DataFrames.missing)
    else
        return DataFrames.DataFrame(wavfile = wf, start = starts, stop = stops, waveform = waveforms)
    end
end


function Plots.plot(click::ClickDetector.ClickPointer)
    plot(ClickDetector.get_data(click)...)
end

function click_spectrum(click::ClickDetector.ClickPointer,
        background_width, n, noverlap=div(n, 2); args...)
    click_samples = ClickDetector.get_data(click)[2]
    background_samples = ClickDetector.get_background(click, background_width)
    pg = welch_pgram(click_samples, n, noverlap; args...)
    pg_back = welch_pgram(background_samples, n, noverlap; args...)
    return(pg, pg_back)
end

function plot_spectrum(click::ClickDetector.ClickPointer, background_width, n,
        noverlap=div(n, 2); args...)
    pg, pg_back = click_spectrum(click, background_width, n, noverlap; args...)
    plot(pg.freq / 1e3, 10log10.(pg.power), xlim=[0, 50], label="Click",
        xlab="Frequency (kHz)", ylab="PSD (dB)")
    plot!(pg_back.freq / 1e3, 10log10.(pg_back.power), xlim=[0, 90], label="Background")
end

end # module

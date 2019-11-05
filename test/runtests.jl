using ClickDetection
using SampledSignals
using DSP
using Statistics
using Dates
using FileIO
using Test

@testset "Intervals" begin
    intervals1 = [1..2, 3..4]
    intervals2 = [1..2, 1.5..4]
    @test length(merge_intervals(intervals1)) == 2
    @test length(merge_intervals(intervals2)) == 1
end

@testset "Detection" begin
    wavfile = joinpath(@__DIR__, "data", "MARS_20180201_131212_clip.wav")
    audio = loadstreaming(wavfile)
    x = read(audio)
    close(audio)

    # preprocess raw audio before click detection
    x = convert.(Float32, x)
    fs = samplerate(x)
    butter = digitalfilter(Highpass(5e3, fs=fs), Butterworth(3))
    x = filtfilt(butter, x)

    # Actual click detection
    thresh = 1e-5
    window = 1e-3
    clicks1 = detect_clicks(x, thresh, window)
    @test series(clicks1[1]) === series(clicks1[2])
    c1 = first(clicks1)
    c2 = Click(c1)
    @test length(samples(c1)) == length(samples(c2))
    @test length(times(c1)) == length(times(c2))
    @test all(samples(c1) .== samples(c2))
    @test all(times(c1) .== times(c2))
    @test left(c1) == left(c2)
    @test right(c1) == right(c2)

    for c in [c1, c2]
        t = times(c)
        @test length(samples(c)) == length(t)
        @test left(c) == first(t)
        @test right(c) == last(t)
    end

    clicks2 = detect_clicks(Array(x), samplerate(x), thresh, window)
    @test series(first(clicks1)) !== series(first(clicks2))
    @test length(clicks1) == length(clicks2)
    @test all(samples(first(clicks1)) .== samples(first(clicks2)))
    @test interval(first(clicks1)) == interval(first(clicks2))
    @test left(first(clicks1)) == left(first(clicks2))
    @test right(first(clicks1)) == right(first(clicks2))
    @test samplerate(first(clicks1)) == samplerate(first(clicks2))
    @test times(first(clicks1)) == times(first(clicks2))

    df = as_dataframe(clicks1, DateTime(2019,1,1,0,0,0), "fakefile.wav")
    @test size(df,1) == length(clicks1)
end

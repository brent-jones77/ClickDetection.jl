using ClickDetection
using SampledSignals
using Statistics
using Test

@testset "ClickDetection.jl" begin
    wavfile = joinpath(@__DIR__, "data", "MARS_20180401_001914.wav")
    x = read_as_axisarray(wavfile)

    noise_pct = 0.6
    noise_factor = 50
    window = 1e-3
    x_teager = teager(x)
    noise_floor = quantile(x_teager, noise_pct)
    thresh = noise_floor * noise_factor

    clicks = detect_clicks(x, thresh, window)
    println(length(clicks))
end

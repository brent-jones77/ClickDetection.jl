module ClickDetection

using AxisArrays
using IntervalSets
using FileIO
using LibSndFile
using SampledSignals
using Statistics
using DataFrames
using Dates

include("detector.jl")
include("io.jl")

export read_as_axisarray,
    AbstractClick,
    ClickPointer,
    Click,
    left,
    right,
    series,
    interval,
    samplerate,
    samples,
    times,
    teager,
    merge_intervals,
    detect_clicks,
    peak_lead_ratio,
    as_dataframe

end # module

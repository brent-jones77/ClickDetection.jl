module ClickDetection

using AxisArrays
using IntervalSets
# using DSP
using FileIO
using LibSndFile
using SampledSignals
using Statistics

include("detector.jl")
include("io.jl")

export read_as_axisarray,
    ClickPointer,
    left,
    right,
    teager,
    merge_intervals,
    detect_clicks,
    get_data,
    peak_lead_ratio

end # module

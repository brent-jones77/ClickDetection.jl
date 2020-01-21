abstract type AbstractClick{T,N} end

struct ClickPointer{T,N} <:  AbstractClick{T,N}
    series::AbstractArray{T,N}
    interval::ClosedInterval
    threshold::Real
end

series(c::ClickPointer) = c.series
interval(c::ClickPointer) = c.interval
left(c::ClickPointer) = interval(c).left
right(c::ClickPointer) = interval(c).right
SampledSignals.samplerate(c::ClickPointer) = samplerate(series(c))
samples(c::ClickPointer) = Array(series(c)[left(c)*s..right(c)*s])
times(c::ClickPointer) = left(c):(1/samplerate(c)):right(c)

struct Click{T,N} <: AbstractClick{T,N}
    samples::AbstractArray{T,N}
    samplerate::Real
    left::Real
    right::Real
end
Click(c::ClickPointer) = Click(samples(c), samplerate(c), left(c),  right(c))

left(c::Click) = c.left
right(c::Click) = c.right
SampledSignals.samplerate(c::Click) = c.samplerate
samples(c::Click) = c.samples
times(c::Click) = left(c):(1/samplerate(c)):right(c)

"""
    teager(x::AbstractArray)

Calculate the Teager energy operator of a sampled signal `x`, in the form of a
`Vector` or n x 1 `Array`.  The Teager energy at sample \$x_i\$ is defined as

```math
x_i^2 - x_{i-1} x_{i+1},
```

so `x` must be at least 3 samples long.
"""
function teager(x::AbstractArray)
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
    merge_intervals(ints::Vector{I}, sort_first=true) where {I <: ClosedInterval}

Combine intervals which overlap.  Takes a vector `ints` of `ClosedIntervals`. If
`sort_first=true` (the default), `ints` is sorted by the left boundary of each
interval.  The only time `sort_first` should be `false` is in the case where `ints`
is already in ascending order.
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

"""
    detect_clicks(x::AbstractArray{T,N}, samplerate, thresh, span)
    detect_clicks(x::AbstractSampleBuf{T,N}, thresh, span)

Detect clicks in a sampled signal using the Teager energy operator. `x` can be an
`AbstractArray` or `AbstractSampleBuf` (from SampledSignals.jl).  If it is the
former, a samplerate must be supplied as well.

`thresh` is the detection threshold for the Teager energy.  `span` is the time
before and after each sample above the Teager threshold to include in the click.
The duration of the shortest click will therefore be `2 * thresh`.

Returns a `Vector{ClickPointer}`.
"""
function detect_clicks(x::AbstractArray{T,N}, samplerate, thresh, span)  where {T, N}
    y = SampleBuf(Array(x), samplerate)
    return detect_clicks(y, thresh, span)
end


function detect_clicks(x::AbstractSampleBuf{T,N}, thresh, span) where {T, N}
    nchannels(x) == 1 && ArgumentError("The SampleBuff cannot have more than one channel.")
    x_teager = teager(x)
    high_teager_idx = findall(x -> x > thresh, x_teager)
    times = domain(x)
    t0, t1 = times[1], times[end]
    intervals = [max(t0, times[i]-span)..min(times[i]+span, t1) for i in high_teager_idx]
    intervals = merge_intervals(intervals)
    return [ClickPointer(x, i, thresh) for i in intervals]
end
#
# function peak_lead_ratio(click::ClickPointer, nsamples_lead)
#     t, x = get_data(click)
#     x² = abs2.(x)
#     return mean(x²[nsamples_lead:end]) / mean(x²[1:nsamples_lead])
# end

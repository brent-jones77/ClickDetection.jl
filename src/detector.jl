
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

"""
Detect clicks in a
"""
function detect_clicks(x::AbstractVector, times::AbstractVector, thresh, span)
    x_teager = teager(x)
    high_teager_idx = findall(x -> x > thresh, x_teager)
    intervals = [max(t0, times[i]-span)..min(times[i]+span, t1) for i in high_teager_idx]
    intervals = merge_intervals(intervals)
    return [ClickPointer(x, i, thresh) for i in intervals]
end

function detect_clicks(x::AxisArray, thresh, span)
    times = axisvalues(x[Axis{:time}])[1]
    t0, t1 = times[1], times[end]
    detect_clicks(x, times, thresh, span)
end

function detect_clicks(x::SampleBuf, thresh, span)
    nchannels(x) == 1 && ArgumentError("The SampleBuff cannot have more than one channel.")
    times = domain(x)
    t0, t1 = times[1], times[end]
    detect_clicks(x, times, thresh, span)
end


function get_data(click::ClickPointer)
    x = click.series[click.interval]
    t = axisvalues(x)[1]
    return t, x
end


function peak_lead_ratio(click::ClickPointer, nsamples_lead)
    t, x = get_data(click)
    x² = abs2.(x)
    return mean(x²[nsamples_lead:end]) / mean(x²[1:nsamples_lead])
end

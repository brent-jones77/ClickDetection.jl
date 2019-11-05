
"""
    as_dataframe(clicks::Vector{T}, starttime::DateTime, wavfile) where {T<:AbstractClick}
    as_dataframe(clicks::Vector{T}, starttime::DateTime) where {T<:AbstractClick}
    as_dataframe(clicks::Vector{T}) where {T<:AbstractClick}

Converts a vector of `ClickPointer` or `Click` structures to tabular form as a
`DataFrame`. If optional argument `starttime` is provided, the `DataFrame`
the clicks' start and end offsets will be converted to `DateTime`s relative
to `starttime` (i.e., the time of the first sample in the analyzed file).  The
source file for the clicks can be included via th optional argument `wavfile`
(either a `String` or `Vector{String}`).
"""
function as_dataframe(clicks::Vector{T}, starttime::DateTime) where {T<:AbstractClick}
    df = DataFrame(
        (start=starttime+Dates.Millisecond(round(left(c) * 1e3)),
        stop=Dates.Millisecond(round(right(c) * 1e3)),
        waveform=samples(c)) for c in clicks)
    return df
end

function as_dataframe(clicks::Vector{T}) where {T<:AbstractClick}
    df = DataFrame((start=left(c), stop=right(c), waveform=samples(c)) for c in clicks)
    return df
end

function as_dataframe(clicks::Vector{T}, starttime::DateTime, wavfile) where {T<:AbstractClick}
    df = as_dataframe(clicks, starttime)
    df[!, :wavfile] .= wavfile
    return df
end

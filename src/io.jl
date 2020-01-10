
"""
    as_dataframe(clicks::Vector{T}, starttime::DateTime, wavfile) where {T<:AbstractClick}
    as_dataframe(clicks::Vector{T}, starttime::DateTime) where {T<:AbstractClick}
    as_dataframe(clicks::Vector{T}) where {T<:AbstractClick}

Converts a vector of `ClickPointer` or `Click` structures to tabular form as a
`DataFrame`. If optional argument `starttime` is provided, the returned
`DataFrame` will also include `start_datetime` and `stop_datetime` columns,
with the clicks' start and end times relative to the start of the sound file
converted to `DateTime`s relative to `starttime` (i.e., the time of the first
sample in the analyzed file).  The source file for the clicks can be included
via th optional argument `wavfile` (either a `String` or `Vector{String}`).
"""
function as_dataframe(clicks::Vector{T}, starttime::DateTime) where {T<:AbstractClick}
    df = DataFrame(
        (start=left(c),
         stop=right(c),
         start_datetime=starttime + Dates.Millisecond(round(left(c) * 1e3)),
         stop_datetime=starttime + Dates.Millisecond(round(right(c) * 1e3)),
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



function as_dataframe(clicks::Vector{T}, starttime, wavfile) where {T<:AbstractClick}
    df = DataFrame(
        (wavfile=wavfile,
        start=starttime+Dates.Millisecond(round(left(c) * 1e3)),
        stop=Dates.Millisecond(round(right(c) * 1e3)),
        waveform=samples(c)) for c in clicks)
    return df
end

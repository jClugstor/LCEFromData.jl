struct LCEProblem
    timeseries::AbstractArray
    embedded_data::AbstractArray
end

struct LCESolution

end

struct EmbeddedData
    timeseries::AbstractArray
    embedded_data::AbstractArray
    dim::Number
    tau::Number

end
struct LCEProblem
    timeseries::Array
    embedded_data::EmbeddedData
end

struct LCESolution
    LCEs::Array
end

struct EmbeddedData
    timeseries::AbstractArray
    embedded_data::AbstractArray
    dim::Number
    tau::Number
end
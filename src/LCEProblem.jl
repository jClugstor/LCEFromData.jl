struct LCEProblem
    timeseries::Array
    timestep::Number
    embedded_data::EmbeddedData
end

struct EmbeddedData
    timeseries::Dataset
    embedded_data::Dataset
    dim::Number
    tau::Number
end

struct LCESolution
    maxLCE::Number
    LCEs::Array
end
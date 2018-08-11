module NoisyNeuron

using Parameters, Random

@with_kw struct TimeAxis
  dt::Float64
  N::Int64
end

include("point_models.jl")

export TimeAxis

end
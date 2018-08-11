module NoisyNeuron

using Parameters, Random

export TimeAxis, PointNeuron, point_white_vtime

@with_kw struct TimeAxis
  dt::Float64
  N::Int64
end

@with_kw struct PointNeuron
  E_L::Float64 # leak reversal potential, mV
  τ_v::Float64 # membrane time constant, ms	
  σ_s::Float64 # noise strength, mV
end

function point_white_vtime(args::PointNeuron, T::TimeAxis; seed = 0)
	if seed != 0
	  Random.seed!(seed)
	end
  @unpack E_L, τ_v, σ_s = args
  @unpack dt, N = T

	V = zeros(N)
  for i=2:N
  	V[i] = V[i-1]*(1-dt/τ_v)+σ_s*sqrt(2*dt/τ_v)*randn()
  end
  t = range(0,step=dt,length=N)
  return t, V .+ E_L
end

end
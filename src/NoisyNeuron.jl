module NoisyNeuron

using Parameters, Random

export TimeAxis, PointNeuron, point_white_vtime, point_white_test

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

function point_white_test()
  T = TimeAxis(dt = 0.1, N = 1000)
  example = PointNeuron(E_L = -60.0,τ_v = 10.0, σ_s = 1.0)
  @time t, V = point_white_vtime(example, T)
  return t, V
end

end
export PointNeuron, WhiteSteady, WhiteCmod, WhiteVmod, point_white_steady_test, point_white_cmod_test, point_white_vmod_test

@with_kw struct PointNeuron
  E_L::Float64 # leak reversal potential, mV
  τ_v::Float64 # membrane time constant, ms	
end

abstract type DriveParams end

@with_kw struct WhiteSteady <: DriveParams
  σ_s::Float64
  E_s::Float64
end

function white_drive(I_args::WhiteSteady, Δ, t)
  @unpack σ_s, E_s = I_args
  E_s*Δ+σ_s*sqrt(2*Δ)*randn()
end

@with_kw struct WhiteCmod <: DriveParams
  σ_s::Float64
  E_s::Float64
  ϵ_c::Float64
  Ω::Float64
end

function white_drive(I_args::WhiteCmod, Δ, t)
  @unpack σ_s, E_s, ϵ_c, Ω = I_args
  E_s*Δ+σ_s*sqrt(2*Δ)*randn()+ϵ_c*Δ*sin(Ω*t)
end

@with_kw struct WhiteVmod <: DriveParams
  σ_s::Float64
  E_s::Float64
  ϵ_v::Float64
  Ω::Float64
end

function white_drive(I_args::WhiteVmod, Δ, t)
  @unpack σ_s, E_s, ϵ_v, Ω = I_args
  E_s*Δ+σ_s*(1+ϵ_v*sin(Ω*t))*sqrt(2*Δ)*randn()
end

function point_white_vtime(N_args::PointNeuron, I_args::DriveParams, T::TimeAxis; seed = 0)
	if seed != 0
	  Random.seed!(seed)
	end
  @unpack E_L, τ_v = N_args
  @unpack dt, N = T

	V = zeros(N)
  for i=2:N
  	V[i] = V[i-1]*(1-dt/τ_v)+white_drive(I_args,dt/τ_v,(i-1)*dt)
  end
  t = range(0,step=dt,length=N)
  return t, V .+ E_L
end

function point_white_steady_test()
  T = TimeAxis(dt = 0.1, N = 1000)
  example = PointNeuron(E_L = -60.0,τ_v = 10.0)
  drive = WhiteSteady(σ_s = 1.0, E_s = 5.0)
  @time t, V = point_white_vtime(example, drive, T; seed = 1000)
  return t, V
end

function point_white_cmod_test()
  T = TimeAxis(dt = 0.1, N = 1000)
  example = PointNeuron(E_L = -60.0,τ_v = 10.0)
  drive = WhiteCmod(σ_s = 1.0, E_s = 5.0, ϵ_c = 5.0, Ω = 0.1*pi)
  @time t, V = point_white_vtime(example, drive, T; seed = 1000)
  return t, V
end

function point_white_vmod_test()
  T = TimeAxis(dt = 0.1, N = 1000)
  example = PointNeuron(E_L = -60.0,τ_v = 10.0)
  drive = WhiteVmod(σ_s = 1.0, E_s = 5.0, ϵ_v = 5.0, Ω = 0.1*pi)
  @time t, V = point_white_vtime(example, drive, T; seed = 1000)
  return t, V
end
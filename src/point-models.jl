export PointNeuron, WhiteSteady, WhiteCmod, WhiteVmod

@with_kw struct PointNeuron
  τ_v::Float64 # membrane time constant, ms	
end

abstract type DriveParams end

@with_kw struct WhiteSteady <: DriveParams
  μ::Float64
  σ_w::Float64
end

function white_drive(I_args::WhiteSteady, dt, t)
  @unpack μ, σ_w = I_args
  μ*dt+σ_w*sqrt(2*dt)*randn()
end

@with_kw struct WhiteCmod <: DriveParams
  μ::Float64
  σ_w::Float64
  ϵ_c::Float64
  Ω::Float64
end

function white_drive(I_args::WhiteCmod, dt, t)
  @unpack μ,σ_w, ϵ_c, Ω = I_args
  μ*dt+σ_w*sqrt(2*dt)*randn()+ϵ_c*dt*sin(Ω*t)
end

@with_kw struct WhiteVmod <: DriveParams
  μ::Float64
  σ_w::Float64
  ϵ_v::Float64
  Ω::Float64
end

function white_drive(I_args::WhiteVmod, dt, t)
  @unpack μ, σ_w, ϵ_v, Ω = I_args
  μ*dt+σ_w*(1+ϵ_v*sin(Ω*t))*sqrt(2*dt)*randn()
end

function point_white_vtime(N_args::PointNeuron, I_args::DriveParams, T::TimeAxis; seed = 0)
	if seed != 0
	  Random.seed!(seed)
	end
  @unpack τ_v = N_args
  @unpack dt, N = T

	V = zeros(N)
  for i=2:N
  	V[i] = V[i-1]*(1-dt)+white_drive(I_args,dt,(i-1)*dt)
  end
  t = range(0,step=dt,length=N)
  return t.*τ_v, V
end
export point_white_vtime

function point_white_steady_test()
  T = TimeAxis(dt = 0.05, N = 1000)
  example = PointNeuron(τ_v = 10.0)
  drive = WhiteSteady(μ=5.0, σ_w = 1.0)
  @time t, V = point_white_vtime(example, drive, T; seed = 1000)
  return t, V
end
export point_white_steady_test

function point_white_cmod_test()
  T = TimeAxis(dt = 0.05, N = 1000)
  example = PointNeuron(τ_v = 10.0)
  drive = WhiteCmod(μ=5.0, σ_w = 1.0, ϵ_c = 5.0, Ω = 0.1*pi)
  @time t, V = point_white_vtime(example, drive, T; seed = 1000)
  return t, V
end
export point_white_cmod_test

function point_white_vmod_test()
  T = TimeAxis(dt = 0.05, N = 1000)
  example = PointNeuron(τ_v = 10.0)
  drive = WhiteVmod(μ=5.0, σ_w = 1.0, ϵ_v = 5.0, Ω = 0.1*pi)
  @time t, V = point_white_vtime(example, drive, T; seed = 1000)
  return t, V
end
export point_white_vmod_test
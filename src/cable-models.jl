export Neurite, WhiteSteadyDist, ColourSteadyDist, sealed_white_var, sealed_white_steady_var_test

@with_kw struct Neurite
  τ_v::Float64
  λ::Float64
  M::Int64 # number of spatial steps
  dx::Float64 # spatial step size
end

abstract type DriveParams end

@with_kw struct WhiteSteadyDist <: DriveParams
  μ::Float64
  σ_w::Float64
end

@with_kw struct ColourSteadyDist <: DriveParams
  μ::Float64
  σ_s::Float64
  τ_s::Float64
end

function white_drive(I_args::WhiteSteadyDist, Δt)
  @unpack σ_w = I_args
  σ_w*sqrt(2*Δt)*randn()
end

function sealed_white_steady_var(N_args::Neurite, I_args::DriveParams, T::TimeAxis; seed=0)
  if seed != 0
    Random.seed!(seed)
  end
  @unpack τ_v,λ,M,dx = N_args
  @unpack μ,σ_w = I_args
  @unpack dt, N = T

  V = zeros(M)
  dV = zeros(M+1)
  ζ = zeros(M)
  V2_sum = zeros(M)
  for i=1:N
    randn!(ζ)
    V[1] = (1-dt)*V[1]+λ^2*(dV[2]-dV[1])dt/dx+2*σ_w*sqrt(λ*dt/dx)*ζ[1]
    V2_sum[1] += V[1]^2
    @simd for j=2:M
      @inbounds V[j] = (1-dt)*V[j]+λ^2*(dV[j+1]-dV[j])*dt/dx+2*σ_w*sqrt(λ*dt/dx)*ζ[j]
      @inbounds dV[j] = (V[j]-V[j-1])/dx
      @inbounds V2_sum[j] += V[j]^2
    end
  end
  return V2_sum/N
end

function sealed_white_steady_var_test()
  T = TimeAxis(dt = 1/500, N = 10000)
  example = Neurite(τ_v = 10.0,λ=200,M=50,dx=20)
  drive = WhiteSteadyDist(μ=5.0,σ_w = 1.0)
  @time Var = sealed_white_steady_var(example, drive, T; seed=1000)
  return Var
end
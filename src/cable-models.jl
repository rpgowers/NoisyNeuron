export Neurite, WhiteSteadyDist, ColourSteadyDist

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
export sealed_white_steady_var

function sealed_colour_steady_var(N_args::Neurite, I_args::DriveParams, T::TimeAxis; seed=0)
  if seed != 0
    Random.seed!(seed)
  end
  @unpack τ_v,λ,M,dx = N_args
  @unpack μ,σ_s,τ_s = I_args
  α_s = τ_s/τ_v
  @unpack dt, N = T
  V2 = zeros(M) # new value of V
  V1 = zeros(M) # old value of V
  S = zeros(M) # synaptic drive term
  dV = zeros(M+1)
  ζ = zeros(M)
  V2_sum = zeros(M)
  dV2_sum = zeros(M)
  for i=1:N
    randn!(ζ)
    V2[1] = V1[1]*(1-dt)+λ^2*(dt/dx)*(dV[2]-dV[1])+S[1]*dt
    S[1] = (1-dt/α_s)*S[1]+2*σ_s*ζ[1]*sqrt(λ*dt/(dx*α_s))
    V2_sum[1] += V2[1]^2
    dV2_sum[1] += ((V2[1]-V1[1])/dt)^2
    V1[1] = V2[1]
    @simd for j=2:M
      @inbounds V2[j] = V1[j]*(1-dt)+λ^2*(dt/dx)*(dV[j+1]-dV[j])+S[j]*dt
      @inbounds S[j] = (1-dt/α_s)*S[j]+2*σ_s*ζ[j]*sqrt(λ*dt/(dx*α_s))
      @inbounds dV[j] = (V2[j]-V2[j-1])/dx
      @inbounds V2_sum[j] += V2[j]^2
      @inbounds dV2_sum[j] += ((V2[j]-V1[j])/dt)^2
      @inbounds V1[j] = V2[j]
    end;
  end
  return V2_sum/N, dV2_sum/N
end
export sealed_colour_steady_var

function sealed_white_steady_var_test()
  T = TimeAxis(dt = 1/500, N = 10000)
  example = Neurite(τ_v = 10.0,λ=200,M=50,dx=20)
  drive = WhiteSteadyDist(μ=5.0,σ_w = 1.0)
  @time Var = sealed_white_steady_var(example, drive, T; seed=1000)
  return Var
end
export sealed_white_steady_var_test

function sealed_colour_steady_var_test()
  T = TimeAxis(dt = 1/500, N = 10000)
  example = Neurite(τ_v = 10.0,λ=200,M=50,dx=20)
  drive = ColourSteadyDist(μ=5.0,σ_s = 1.0,τ_s=5.0)
  @time Var, dVar = sealed_colour_steady_var(example, drive, T; seed=1000)
  return Var, dVar
end
export sealed_colour_steady_var_test
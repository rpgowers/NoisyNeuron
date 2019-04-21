export Neurite, WhiteSteadyDist, ColourSteadyDist

@with_kw struct Neurite
  τ_v::Float64
  λ::Float64
  M::Int64 # number of spatial steps
  dx::Float64 # spatial step size
end

function x_axis(args_N::Neurite)
  @unpack M,dx = args_N
  x = 0.5*dx:dx:M*dx
end
export x_axis

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

function cfunc(x,ζ,L,λ)
  a = sqrt(1+ζ)
  return (2/a)*cosh((L-x)*a/λ)*cosh(x*a/λ)/sinh(L*a/λ)
end

function sealed_steady_var(N_args::Neurite, I_args::WhiteSteadyDist, x)
  @unpack λ,M,dx = N_args
  @unpack σ_w = I_args
  L = M*dx
  return σ_w^2 .* cfunc.(x,0,L,λ)
end
export sealed_steady_var

function sealed_steady_var(N_args::Neurite, I_args::ColourSteadyDist, x)
  @unpack λ,M,dx,τ_v = N_args
  @unpack σ_s,τ_s = I_args
  α_s = τ_s/τ_v
  L = M*dx
  return σ_s^2*α_s.*(cfunc.(x,0,L,λ).-cfunc.(x,1/α_s,L,λ))
end
export sealed_steady_var

function sealed_steady_dvar(N_args::Neurite, I_args::ColourSteadyDist, x)
  @unpack λ,M,dx,τ_v = N_args
  @unpack σ_s,τ_s = I_args
  α_s = τ_s/τ_v
  L = M*dx
  return σ_s^2 .* cfunc.(x,1/α_s,L,λ)./α_s
end
export sealed_colour_steady_dvar

function sealed_steady_var(N_args::Neurite, I_args::WhiteSteadyDist, T::TimeAxis; seed=0)
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
export sealed_steady_var

function sealed_steady_var(N_args::Neurite, I_args::ColourSteadyDist, T::TimeAxis; seed=0)
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
export sealed_steady_var

function sealed_white_steady_var_test()
  T = TimeAxis(dt = 1/500, N = 5000000)
  example = Neurite(τ_v = 10.0,λ=200,M=50,dx=20)
  drive = WhiteSteadyDist(μ=5.0,σ_w = 1.5)
  @time Var_sim = sealed_steady_var(example, drive, T; seed=1000)
  x = x_axis(example)
  Var_th = sealed_steady_var(example,drive,x)
  return Var_sim,Var_th,x
end
export sealed_white_steady_var_test

function sealed_colour_steady_var_test()
  T = TimeAxis(dt = 1/500, N = 5000000)
  example = Neurite(τ_v = 10.0,λ=200,M=50,dx=20)
  drive = ColourSteadyDist(μ=5.0,σ_s = 1.5,τ_s=5.0)
  @time Var_sim, dVar_sim = sealed_steady_var(example, drive, T; seed=1000)
  x = x_axis(example)
  Var_th = sealed_steady_var(example,drive,x)
  dVar_th = sealed_steady_dvar(example,drive,x)
  return Var_sim, dVar_sim, Var_th, dVar_th, x
end
export sealed_colour_steady_var_test
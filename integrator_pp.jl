using DifferentialEquations, Distributed
using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit
norminf(x) =  norm(x, Inf)

function pp!(du, u, p, t)
  @unpack beta1, beta2, d, h, c, zeta = p
  x, y, z = u
  sig1 = x + beta1
  sig2 = x + beta2

 du[1] = x * (1 - x - y/sig1 - z/sig2)
 du[2] = zeta * y * (x/sig1 - c)
 du[3] = zeta * z * (x/sig2 - d - h*z)
 du
end

par = (beta1 = 0.2, beta2 = 0.25, d = 0.17, h = 1.3, c = 0.2, zeta = 0.11)
u0 = [0.3, 0.4, 0.1]
tspan = (0.0, 800.0)
prob = ODEProblem(pp!, u0, tspan, par)
h = 0.01
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot(sol, plotdensity=5000)


@everywhere function prob_func(prob,i,repeat)
  remake(prob, u0=rand()*prob.u0)
end

traj_no = 8
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sim = solve(ensemble_prob, Tsit5(), EnsembleDistributed(), trajectories=8)
plot(sim, vars=(1,2,3), plotdensity=1000)
plot!([u0[1]], [u0[2]], [u0[3]])
for i in 1:traj_no
  x0 = sim[i].prob.u0
  scatter!([x0[1]], [x0[2]], [x0[3]], markersize=3)
end
title!("Predator-Prey model: h = 1.3")
savefig("pp_h1.3.png")

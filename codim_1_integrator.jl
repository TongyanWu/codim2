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

par = (beta1 = 0.2, beta2 = 0.25, d = 0.17, h = 1.02, c = 0.2, zeta = 0.11)
u0 = [0.15, 0.4, 0.25]
tspan = (0.0, 8000.0)
prob = ODEProblem(pp!, u0, tspan, par)
h = 0.01
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot(sol, vars =(1,2,3), plotdensity=5000, labels="")
scatter!([u0[1]], [u0[2]], [u0[3]], c=:lightblue, labels="(0.3, 0.4, 0.1)")
u0 = [0.15, 0.4, 0.02]
tspan = (0.0, 2000.0)
prob = ODEProblem(pp!, u0, tspan, par)
h = 0.01
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot!(sol, vars =(1,2,3), plotdensity=5000, c=:green, labels="")
scatter!([u0[1]], [u0[2]], [u0[3]], c=:green,labels="(0.15, 0.4, 0.02)")
u0 = [0.3363, 2.659e-22, 0.38912]
scatter!([u0[1]], [u0[2]], [u0[3]], c=:orange,labels="(0.336, 2.66e-22, 0.389)")
#savefig("pp_outside_109.png")
# third initial value
#tspan = (0.0, 200.0)
#prob = ODEProblem(pp!, u0, tspan, par)
#h = 0.01
#sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
#plot(sol, vars =(1,2,3), plotdensity=5000, c=:orange, labels="")
#scatter!([u0[1]], [u0[2]], [u0[3]], c=:orange,labels="(0.336, 2.66e-22, 0.389)")
#title!("System (1.1): h = 1.04")
#savefig("pp_inside_109.png")

# fourth init
u0 = [0.5, 2.659e-22, 0.5]
tspan = (0.0, 800.0)
prob = ODEProblem(pp!, u0, tspan, par)
h = 0.01
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot(sol, vars =(1,2,3), plotdensity=5000, c=:red, labels="")
scatter!([u0[1]], [u0[2]], [u0[3]], c=:red,labels="(0.5, 2.66e-22, 0.5)")

#title!("System (1.1): h = 1.04")

@everywhere function prob_func(prob,i,repeat)
  remake(prob, u0=rand()*prob.u0)
end

traj_no = 15
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sim = solve(ensemble_prob, Tsit5(), EnsembleDistributed(), trajectories=8)
plot(sim, vars=(1,2,3), plotdensity=1000)

for i in 1:traj_no
  x0 = sim[i].prob.u0
  scatter!([x0[1]], [x0[2]], [x0[3]], markersize=3)
end
scatter!([u0[1]], [u0[2]], [u0[3]])
title!("Predator-Prey model: h = 1.037")
savefig("pp_h1037.png")

## h = 1.0
par = (beta1 = 0.2, beta2 = 0.25, d = 0.17, h = 1.0, c = 0.2, zeta = 0.11)
tspan = (0.0, 130.0)
prob = ODEProblem(pp!, u0, tspan, par)
traj_no = 8
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sim = solve(ensemble_prob, Tsit5(), EnsembleDistributed(), trajectories=8)
plot(sim, vars=(1,2,3), plotdensity=1000)
plot!([u0[1]], [u0[2]], [u0[3]])

for i in 1:traj_no
  x0 = sim[i].prob.u0
  scatter!([x0[1]], [x0[2]], [x0[3]], markersize=3)
end
title!("Predator-Prey model: h = 1.0")
savefig("pp_h1.png")
### one solution
u0 = [0.26, 0.01, 0.15]
par = (beta1 = 0.2, beta2 = 0.25, d = 0.17, h = 0.1, c = 0.2, zeta = 0.11)
tspan = (0.0, 4030.0)
prob = ODEProblem(pp!, u0, tspan, par)
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot(sol, vars=(1,2,3), plotdensity=3000)
scatter!([u0[1]], [u0[2]], [u0[3]])

#### h1 only one

par = (beta1 = 0.2, beta2 = 0.25, d = 0.17, h = 0.9, c = 0.2, zeta = 0.11)
u0 = [0.3, 0.4, 0.1]
tspan = (0.0, 8000.0)
prob = ODEProblem(pp!, u0, tspan, par)
h = 0.01
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot(sol, vars =(1,2,3), plotdensity=5000)
scatter!([u0[1]], [u0[2]], [u0[3]],labels="initial value 1")

u0 = [0.15, 0.4, 0.02]
tspan = (0.0, 8000.0)
prob = ODEProblem(pp!, u0, tspan, par)
h = 0.01
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot!(sol, vars =(1,2,3), plotdensity=5000, c=:green)
scatter!([u0[1]], [u0[2]], [u0[3]],labels="initial value 2")
title!("h=0.9")
savefig("h09.png")



par = (beta1 = 0.2, beta2 = 0.25, d = 0.17, h = 1.09, c = 0.2, zeta = 0.11)
#1st initial value
u0 = [0.3, 0.4, 0.1]
tspan = (0.0, 2000.0)
prob = ODEProblem(pp!, u0, tspan, par)
h = 0.01
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot(sol, vars =(1,2,3), plotdensity=5000, labels="")
scatter!([u0[1]], [u0[2]], [u0[3]], c=:lightblue, labels="(0.3, 0.4, 0.1)")
#2nd
u0 = [0.15, 0.4, 0.02]
tspan = (0.0, 2000.0)
prob = ODEProblem(pp!, u0, tspan, par)
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot!(sol, vars =(1,2,3), plotdensity=5000, c=:green, labels="")
scatter!([u0[1]], [u0[2]], [u0[3]], c=:green,labels="(0.15, 0.4, 0.02)")

# 3rd (y=2.66e-22 or y= 0)
u0 = [0.3363, 2.66e-22 , 0.38912]
scatter([u0[1]], [u0[2]], [u0[3]], c=:orange,labels="(0.336, 2.66e-22 , 0.389)",xlabel="x", ylabel="y", zlabel="z")
# third initial value
tspan = (0.0, 920.0)
prob = ODEProblem(pp!, u0, tspan, par)
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot!(sol, vars =(1,2,3), plotdensity=5000, c=:orange, labels="")
#scatter!([u0[1]], [u0[2]], [u0[3]], c=:orange,labels="(0.336, 2.66e-22, 0.389)")
#title!("System (1.1): h = 1.04")
#savefig("pp_inside_109.png")

# 4th init
u0 = [0.5, 2.66e-22 , 0.5]
tspan = (0.0, 1200.0)
prob = ODEProblem(pp!, u0, tspan, par)
h = 0.01
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot!(sol, vars =(1,2,3), plotdensity=5000, c=:red, labels="")
scatter!([u0[1]], [u0[2]], [u0[3]], c=:red,labels="(0.5, 2.66e-22 , 0.5)")

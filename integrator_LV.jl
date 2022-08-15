using DifferentialEquations
using Plots
using Distributed
### Plotting Chapter 2 Lotka Volterra Phase-Portrait

# ODE Lotka Volterra Classic

function lv!(dx, x, p, t)
 dx[1] = x[1] - x[1].*x[2]
 dx[2] = -x[2] + x[1].*x[2]
end

x0 = [1.0, 2.0]
tspan = (0.0, 500.0)
prob = ODEProblem(lv!, x0, tspan)
h = 0.01
# try also ImplicitEurler() or Euler()
sol = solve(prob,  Midpoint(), dt=h, adaptive=false)
plot(sol, vars=(1,2), plotdensity=500, xlabel="Prey x",
    ylabel="Predator y", labels="(x(t), y(t))")
title!("Lotka-Volterra: Classical")
#savefig("LV.png")

#init = [r for r in 1:10, c in 1:2]
init = [1:0.2:2.8  1:0.2:2.8]
initial_conditions = tuple.(init[:, 1], init[:,2])
@everywhere function prob_func(prob,i,repeat)
  (prob,u0=rand()*prob.u0)
end

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sim = solve(ensemble_prob, Midpoint(), EnsembleDistributed(), trajectories=10)
plot(sim, vars=(1,2), plotdensity=500, xlabel="Prey x",
    ylabel="Predator y", labels="(x(t), y(t))")
title!("Lotka-Volterra: Classical")
#savefig("LV_ens.png")
# ODE Lotka Volterra Logistic Growth
function lv_log!(dx, x, p, t)
 dx[1] = x[1]  .* (1 + x[1] / 2) - x[1].*x[2]
 dx[2] = -x[2] + x[1].*x[2]
end

x0 = [1.0, 2.0]
tspan = (0.0, 55.0)
prob = ODEProblem(lv_log!, x0, tspan)
h = 0.01
# try also ImplicitEurler() or Euler()
sol = solve(prob,  Midpoint(), dt=h, adaptive=false)
plot(sol, vars=(1,2), plotdensity=5000, xlabel="Prey x",
    ylabel="Predator y", labels="(x(t), y(t))")
title!("Lotka-Volterra: Logistic Growth")
#savefig("LV_log.png")

# make ensemble (multiple initial values)
init = [1:0.2:2.8  0.9:-0.05:0.45 ]
initial_conditions = tuple.(init[:, 1], init[:,2])
@everywhere function prob_func(prob,i,repeat)
  remake(prob,u0=rand()*prob.u0)
end

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sim = solve(ensemble_prob, Midpoint(), EnsembleDistributed(), trajectories=10)
plot(sim, vars=(1,2), plotdensity=500, xlabel="Prey x",
    ylabel="Predator y", labels="(x(t), y(t))")
title!("Lotka-Volterra: Logistic Growth")
#savefig("LV_log_ens.png")

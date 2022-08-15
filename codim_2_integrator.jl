### This module needs to have integrator_pp.jl loaded first

# Used for plotting Chapter 5.2 phase portraits, trying to find meaningful
# examples showcasing different behavior in Bautin and Fold-Hopf

par = (beta1 = 0.2, beta2 = 1.5, d = 0.17, h = -0.1, c = 0.2, zeta = 0.11)
#1st initial value
u0 = [0.3, 0.4, 0.1]
tspan = (0.0, 1000.0)
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
#par = (beta1 = 0.2, beta2 = 0.34, d = 0.17, h = 1.0, c = 0.2, zeta = 0.11)
u0 = [0.336, 0, 0.389]
scatter([u0[1]], [u0[2]], [u0[3]], c=:orange,labels="(0.336, 0, 0.389)")
# third initial value
tspan = (0.0, 9000.0)
prob = ODEProblem(pp!, u0, tspan, par)
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot!(sol, vars =(1,2,3), plotdensity=5000, c=:orange, labels="", xlabel="x", ylabel="y", zlabel="z")
#scatter!([u0[1]], [u0[2]], [u0[3]], c=:orange,labels="(0.336, 2.66e-22, 0.389)")
#title!("System (1.1): h = 1.04")
#savefig("pp_inside_109.png")

# 4th init
u0 = [0.5, 0.0, 0.5]
tspan = (0.0, 4000.0)
prob = ODEProblem(pp!, u0, tspan, par)
h = 0.01
sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
plot!(sol, vars =(1,2,3), plotdensity=5000, c=:red, labels="")
scatter!([u0[1]], [u0[2]], [u0[3]], c=:red,labels="(0.5, 0, 0.5)")

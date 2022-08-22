using DifferentialEquations
using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
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

######## Find regions of bautin bifurcation
# 1st build cyclic iterator around bifurcation point
coordinates = []
n = 32
eps = 0.1
for k in 1:n
    z = eps * exp(2*pi*im/n)^k
    push!(coordinates,z)
end
# step 2: add values to bifurcation point (h = 1.472222, beta2 = 0.176)
parameters = []
for k in 1:n
    b2 = 0.17667440497694786 + real(coordinates[k])
    h2 = 1.4722229970610037 + imag(coordinates[k])
    par = (beta1 = 0.2, beta2 = b2, d = 0.17, h = h2, c = 0.2, zeta = 0.11)
    push!(parameters, par)
    #create_plot(par)
end

function create_plot(par, counter, n)
    # IV close to equilibrium
    u0 = [0.370471722674988, 2.522332479229523e-27, 0.3444439591866704]
    st = "$u0"
    tspan = (0.0, 500.0)
    prob = ODEProblem(pp!, u0, tspan, par)
    h = 0.01
    sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
    plot(sol, vars =(1,3), plotdensity=5000, labels="")
    #scatter!([u0[1]], [u0[2]], [u0[3]], c=:lightblue, labels="")
    scatter!([u0[1]], [u0[3]], c=:lightblue,labels="")

    #2nd IV further away
    u1 = [0.42, 0.0, 0.4]
    st = "$u1"
    prob = ODEProblem(pp!, u1, tspan, par)
    sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
    #plot!(sol, vars =(1,2,3), plotdensity=5000, c=:green, labels="")
    plot!(sol, vars =(1,3), plotdensity=5000, c=:green, labels="")
    #scatter!([u1[1]], [u1[2]], [u1[3]], c=:green,labels="")
    scatter!([u1[1]],[u1[3]], c=:green,labels="")

    #3rd IV
    u2 = [1.05, 0.0, 1.15]
    st = "$u2"
    prob = ODEProblem(pp!, u2, tspan, par)
    sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
    #plot!(sol, vars =(1,2,3), plotdensity=5000, c=:orange, labels="")
    plot!(sol, vars =(1,3), plotdensity=5000, c=:orange, labels="")
    #scatter!([u2[1]], [u2[2]], [u2[3]], c=:orange,labels="")
    scatter!([u2[1]], [u2[3]], c=:orange,labels="")

    u3 = [2.7, 0.0, 1.9]
    st = "$u3"
    prob = ODEProblem(pp!, u3, tspan, par)
    sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
    #plot!(sol, vars =(1,2,3), plotdensity=5000, c=:red, labels="")
    plot!(sol, vars =(1,3), plotdensity=5000, c=:red, labels="")
    #scatter!([u3[1]], [u3[2]], [u3[3]], c=:red,labels="")
    scatter!([u3[1]], [u3[3]],c=:red, labels="")

    #u4 = [6.31, 0.0, 6.25]
    #prob = ODEProblem(pp!, u4, tspan, par)
    #sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
    #plot!(sol, vars =(1,3), plotdensity=5000, c=:violet, labels="")
    #scatter!([u4[1]], [u4[2]], [u4[3]],labels="")
    #scatter!([u4[1]], [u4[3]],c=:violet, labels="")
    ## very far out
    u5 = [3.5, 0.0, 3.4]
    prob = ODEProblem(pp!, u5, tspan, par)
    sol = solve(prob,  Tsit5(), dt=h, adaptive=false)
    #plot!(sol, vars =(1,2,3), plotdensity=5000, labels="")
    plot!(sol, vars =(1,3), plotdensity=5000,c=:violet, labels="")

    title!("\$ \\epsilon = $eps,\\; \\theta = (\\frac{2\\pi}{$n})^{$counter} \$")
    #scene = scatter!([u5[1]], [u5[2]], [u5[3]],labels="")
    scene = scatter!([u5[1]], [u5[3]],c=:violet,labels="")
    display(scene)
    savefig("$eps _$n _$counter.png")
    #sleep(0.25)
    ### title
    #h_st = par[4]
    #b2_st = par[2]
    #title!("h=$h_st, b2=$b2_st")
end


#j=14
#par = parameters[j]
#create_plot(par)
#title!("(2pi/16)^$j")
for j in 1:n
    par = parameters[j]
    create_plot(par, j, n)
end


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

pp(z, p) = pp!(similar(z), z, p, 0)
dpp(z,p) = ForwardDiff.jacobian(x -> pp(x,p), z)
jet = BK.getJet(pp, dpp)

par = (beta1 = 0.2, beta2 = 0.25, d = 0.17, h = 0.9, c = 0.2, zeta = 0.11)
u0 = [200.0, 50.0, 5.0]

opts_br = ContinuationPar(pMin = 0.1, pMax = 15.0,
               # parameters to have a smooth result
               ds = 0.04, dsmax = 0.05,
               # this is to detect bifurcation points precisely with bisection
               detectBifurcation = 3,
               # Optional: bisection options for locating bifurcations
               nInversion = 8, maxBisectionSteps = 25, nev = 5)

br, = continuation(pp, dpp, u0, par, (@lens _.h), opts_br;
recordFromSolution = (u, p) -> ( x = u[1], y = u[2], z = u[3]),
tangentAlgo = BorderedPred(),
plot = true, normC = norminf)

scene = plot(br, plotfold=false, markersize=3, legend=:topleft)
savefig("continuation.png")


######

diagram = bifurcationdiagram(jet...,
	# initial point and parameter
	u0, par,
	# specify the continuation parameter
	(@lens _.h),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	3,
	(args...) -> setproperties(opts_br; pMin = 0.1, pMax = 15.0, ds = 0.04, dsmax = 0.05, nInversion = 8, detectBifurcation = 3, dsminBisection =1e-18, maxBisectionSteps=20);
	recordFromSolution = (u, p) -> ( x = u[1], y = u[2], z = u[3]) )

	#######

	# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  maxIter = 10)

opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, pMax = 5., pMin=-5.0,
	maxSteps = 140, newtonOptions = (@set optn_po.tol = 1e-7),
	nev = 3, precisionStability = 1e-8, detectBifurcation = 3, plotEveryStep = 20, saveSolEveryStep=1)


	args_po = (	recordFromSolution = (u, p) -> begin
			xtt = BK.getPeriodicOrbit(p.prob, u, @set par.h = p.p)
			return (max = maximum(xtt[1,:]),
					min = minimum(xtt[1,:]),
					period = getPeriod(p.prob, u, @set par.h = p.p))
		end,
		plotSolution = (u, p; k...) -> begin
			xtt = BK.getPeriodicOrbit(p.prob, u, @set par.h = p.p)
			plot!(xtt.t, xtt[1,:]; label = "x", k...)
			plot!(xtt.t, xtt[2,:]; label = "y", k...)
			plot!(xtt.t, xtt[3,:]; label = "z", k...)
			plot!(br; subplot = 1, putspecialptlegend = false)
			end,
		normC = norminf)

		######

		Mt = 200 # number of time sections
	br_potrap, utrap = continuation(jet...,
	# we want to branch form the 2nd st bif. point
	br, 1, opts_po_cont,
	# we want to use the Trapeze method to locate PO
	PeriodicOrbitTrapProblem(M = Mt);
	# this jacobian is specific to ODEs
	# it is computed using AD and
	# updated inplace
	jacobianPO = :Dense,
	# regular continuation options
	verbosity = 2,	plot = true,
	args_po...)

	scene = plot(br, br_potrap, markersize = 3)
plot!(scene, br_potrap.param, br_potrap.min, label = "")
savefig("diagram.png")
##
plot()
# fetch the saved solutions
for sol in br_potrap.sol[1:2:40]
	# periodic orbit
	po = sol.x
	# get the mesh and trajectory
	traj = BK.getPeriodicOrbit(br_potrap.functional, po, @set par.h = sol.p)
	plot!(traj[1,:], traj[3,:], xlabel = "x", ylabel = "z", label = "")
end
title!("Periodic Orbits")
savefig("periodic_orbits.png")
###### codim 2 bifurcation search


hp_codim2, = continuation(jet[1:2]..., br, 1, (@lens _.beta2),
	ContinuationPar(opts_br, pMin = -0.01, pMax = 1.4,
		ds = -0.001, dsmax = 0.01) ;
	normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detectCodim2Bifurcation = 2,
	# this is required to detect the bifurcations
	d2F = jet[3], d3F = jet[4],
	# tell to start the Hopf problem using eigen elements: compute left eigenvector
	startWithEigen = true,
	# we save the first component for plotting
	recordFromSolution = (u,p; kw...) -> (x = u.u[1] ),
	# we update the Hopf problem at every continuation step
	updateMinAugEveryStep = 1,
	# compute both sides of the initial condition
	bothside = true,
	# use this linear bordered solver, better for ODEs
	bdlinsolver = MatrixBLS(),
	)

	scene = plot(hp_codim2, vars=(:beta2, :x), branchlabel = "HOPF")
	plot!(scene, br, xlims=(0.8,1.3))
	plot(hp_codim2)
	savefig("codim2.png")

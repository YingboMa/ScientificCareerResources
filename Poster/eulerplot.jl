using Plots, OrdinaryDiffEq
pgfplots()
prob = ODEProblem((u,p,t)->-50(u-cos(t)), 0.0, (0, 1.5))
imeuler = solve(prob, ImplicitEuler())
imeuler = solve(prob, ImplicitEuler(), tstops=imeuler.t[1:20:end], adaptive=false, dense=false)
exeuler = solve(prob, Euler(), dt=1.9/50, dense=false)
truesol = solve(prob, Vern9(), abstol=1e-10, reltol=1e-10)
plt1 = plot(imeuler, ylims=(0, 3), lab="Implicit Euler", marker=(:hex, 2, 0.8, Plots.stroke(3, :gray)));
plt2 = plot!(plt1, exeuler, lab="Explicit Euler", marker=(:d, 1, 0.8, Plots.stroke(3, :blue)), line = (:dot, 2));
plt3 = plot!(truesol, line=(:path, 1), lab="True Solution", title="Implicit Euler vs Explicit Euler", size=(330, 270), xlabel=nothing, ylabel=nothing, titlefontsize=8, legendfontsize=3);
savefig("stiffness_euler.tex")
#savefig("stiffness_euler.pdf")

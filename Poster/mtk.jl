using ModelingToolkit, OrdinaryDiffEq
@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)
eqs = [D(x) ~ σ*(y-x),D(y) ~ x*(ρ-z)-y,D(z) ~ x*y - β*z]
@named lorenz1 = ODESystem(eqs)
u0 = [x => 1.0,y => 0.0,z => 0.0]
p  = [σ => 10.0, ρ => 28.0, β => 8/3]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz1,u0,tspan,p)
sol = solve(prob,Tsit5())

@named lorenz2 = ODESystem(eqs)
@variables a; @parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + a*γ]
@named connected = ODESystem(connections,t,[a],[γ],systems=[lorenz1,lorenz2])
u0 = [lorenz1.x => 1.0,lorenz1.y => 0.0,lorenz1.z => 0.0,
      lorenz2.x => 0.0,lorenz2.y => 1.0,lorenz2.z => 0.0,
      a => 2.0]
p  = [lorenz1.σ => 10.0,lorenz1.ρ => 28.0,lorenz1.β => 8/3,
      lorenz2.σ => 10.0,lorenz2.ρ => 28.0,lorenz2.β => 8/3,
      γ => 2.0]
tspan = (0.0,100.0)
prob = ODEProblem(connected,u0,tspan,p)
sol = solve(prob,Rodas4())

function pendulum!(du, u, p, t)
    x, dx, y, dy, T = u
    g, L = p
    du[1] = dx; du[2] = T*x
    du[3] = dy; du[4] = T*y - g
    du[5] = x^2 + y^2 - L^2
end
pendulum_fun! = ODEFunction(pendulum!, mass_matrix=Diagonal([1,1,1,1,0]))
u0 = [1.0, 0, 0, 0, 0]; p = [9.8, 1]; tspan = (0, 10.0)
pendulum_prob = ODEProblem(pendulum_fun!, u0, tspan, p)
traced_sys = modelingtoolkitize(pendulum_prob)
pendulum_sys = structural_simplify(dae_index_lowering(traced_sys))
prob = ODAEProblem(pendulum_sys, Pair[], tspan)
sol = solve(prob, Tsit5())

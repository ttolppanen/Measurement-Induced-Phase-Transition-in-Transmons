using MIPTM, Plots, ParametersModule
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl", "OllisCode/Basis.jl", "OllisCode/Entropy.jl"])


function f()
	L = 6; N = L;
	sp = SystemParameters(L=L, N=N, cap=N)
	state = onesState(sp)
	projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=0.1, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters(sp=sp, U=5.14)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=10, dt=0.1, time=30.0, traj=10)
	prop = propagator(p.bhp.𝐻, p.Ψ₀, p.sdim, p.t.dt)
	pic = plot(abs2.(prop), ylims=(0, 0.1))
	display(pic)
	display(abs2.(prop)[8])
end

f()

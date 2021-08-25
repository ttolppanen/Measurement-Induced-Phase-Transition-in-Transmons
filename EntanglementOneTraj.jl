using Plots, MIPTM, ParametersModule
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function f(prob)
	L = 6; N = L;
	cap = N
	sp = SystemParameters(L=L, N=N, cap=cap)
	state = onesState(sp)
	projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=prob, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters(L=L, N=N, cap=cap, U=-0.14)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=3, dt=0.02, time=10.0, traj=1)
	@time sol = MIPT(p, projectAfterTimeStep = false)
	#@time res = calcMean(sol, Ψ->entanglement_entropy(p.sp.L, p.sp.N, Ψ, 2, cap=cap))
	#popfirst!(res)
	#popfirst!(p.t.times)
	#plot(p.t.times, res, label="p=$(p.pp.p)", legend=:left)
end

f(0.0)
#f(0.06)
#f(0.16)

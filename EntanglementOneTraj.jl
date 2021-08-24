using Plots, MIPTM, ParametersModule
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function f(prob)
	L = 4; N = L;
	cap = N
	state = onesState(L, cap)
	projOp = generateProjectionOperators(L, N, cap)
	sp = SystemParameters(L=L, N=N, Ψ₀=state)
	pp = ProjectionParameters(p=prob, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters(L=L, N=N, cap=cap, U=-0.14)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, cap=cap, sdim=3, dt=0.02, time=10.0, traj=10, useKrylov=false)
	@time sol = MIPT(p, projectAfterTimeStep = true)
	#@time res = calcMean(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, 2, cap))
	#popfirst!(res)
	#popfirst!(p.t.times)
	#plot(p.t.times, res, label="p=$(p.p)", legend=:left)
end

f(0.0)
#f(0.06)
#f(0.16)

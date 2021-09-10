using Plots, MIPTM, ParametersModule
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function f(prob)
	L = 10; N = L;
	cap = 2
	sp = SystemParameters(L=L, N=N, cap=cap)
	state = onesState(sp)
	projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=prob, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters(sp=sp, U=-0.14, Uσ=1.07)
	@time p = Parameters(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=3, dt=0.02, time=30.0, traj=1)
	@time sol = MIPT(p, projectAfterTimeStep = true)
	#@time res = calcMean(sol, Ψ->entanglement_entropy(p.sp.L, p.sp.N, Ψ, 3, cap=cap))
	#popfirst!(res)
	#popfirst!(p.t.times)
	#plot(p.t.times, res, label="p=$(p.pp.p)", legend=:left)
end

@time f(0.06)
#f(0.06)
#f(0.16)

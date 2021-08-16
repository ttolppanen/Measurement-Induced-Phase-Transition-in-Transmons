using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function f(prob)
	L = 8; N = L;
	cap = N - 6
	state = onesState(L, cap)
	measOp = generateProjectionOperators(L, N, cap)
	p = ParametersConstructor(L=L, N=N, cap=cap, sdim=3, measOp=measOp, traj=4, dt=0.02, time=3.0, p=prob, f=1.0, U=0.14, J=1.0, Ψ₀=state)
	@time sol = MIPT(p, projectAfterTimeStep = false)
	#@time res = calcMean(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, 2, cap))
	#popfirst!(res)
	#popfirst!(p.t.times)
	#plot(p.t.times, res, label="p=$(p.p)", legend=:left)
end

f(0.0)
#f(0.06)
#f(0.16)

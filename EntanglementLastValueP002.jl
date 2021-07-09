using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function f(prob)
	L = 8; N = L;
	state = onesState(L)
	measOp = generateProjectionOperators(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=10, measOp=measOp, traj=100, dt=0.02, time=10.0, p=prob, f=1.0, U=0.14, J=1.0, Ψ₀=state)
	@time sol = MIPT(p, projectAfterTimeStep = true, onlyLastValue = true)
	res = calcMean(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, 4))[1]
	display(res)
end

f(0.02)

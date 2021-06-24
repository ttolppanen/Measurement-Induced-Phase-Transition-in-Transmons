using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function f()
	L = 8; N = 4;
	state = zeros(dimension(L, N))
	state[find_index([1, 0, 1, 0, 1, 0, 1, 0])] = 1.;
	measOp = generateProjectionOperators(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=10, measOp=measOp, dt=0.01, time=10.0, p=0.1, f=1.0, U=1.0, J=1.0, Ψ₀=state)
	sol = solveEveryTimeStep(p)
	res = calcMean([sol], Ψ->entanglement_entropy(p.L, p.N, Ψ, 4))
	plot(p.t.times, res)
end

@time f()

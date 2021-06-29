using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function calcWithDifferentProb(L, prob)
	N = L
	state = zeros(dimension(L, N))
	state[find_index(ones(Int64, L))] = 1.;
	measOp = singleSubspaceProjectors(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=10, measOp=measOp, traj=50, dt=0.1, time=10.0, p=prob, f=10.0, U=0.14, J=1.0, Ψ₀=state)
	@time sol = MIPT(p)
	@time res = calcMean(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, Int(L / 2)))
	return res, p
end

function f()
	L = 6
	prob = [0.01, 0.06, 0.1, 0.2, 0.4, 0.8]
	res, param = calcWithDifferentProb(L, 0.0)
	pl = plot(param.t.times, res, label="p = 0.0")
	for p in prob
		res, param = calcWithDifferentProb(L, p)
		plot!(param.t.times, res, label="p = $p")
	end
	return pl
end

f()

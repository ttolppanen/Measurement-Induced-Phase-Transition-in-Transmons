using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function f()
	L = 6; N = L;
	state = zeros(dimension(L, N))
	state[find_index(ones(Int64, L))] = 1.;
	measOp = singleSubspaceProjectors(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=10, measOp=measOp, traj=100, dt=0.1, time=10.0, p=0.50, f=1.0, U=1.0, J=1.0, Ψ₀=state)
	@time sol = MIPT(p)
	@time res = calcMean(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, Int(L / 2)))
	plot(p.t.times, res)
end

f()

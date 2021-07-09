using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function f(prob)
	L = 6; N = 6;
	state = zeros(dimension(L, N))
	state[find_index([1, 1, 1, 1, 1, 1])] = 1.;
	display(dimension(L, N))
	measOp = singleSubspaceProjectors(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=10, measOp=measOp, traj=100, dt=0.02, time=3.0, p=prob, f=1.0, U=0.14, J=1.0, Ψ₀=state)
	@time sol = MIPT(p, projectAfterTimeStep = true)
	@time res = calcMean(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, 3))
	popfirst!(res)
	popfirst!(p.t.times)
	plot!(p.t.times, res, xaxis=:log, label="p=$(p.p)")
end

#f(0.01)
f(0.06)
f(0.16)

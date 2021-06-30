using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function calcWithDifferentProb(p::Parameters, probabilities)
	out = []
	for prob in probabilities
		param = ParametersConstructorWithP(p, prob)
		@time sol = MIPTProjectAfterEveryTimeStep(param)
		@time res = calcMean(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, Int(p.L / 2)))
		push!(out, res)
	end
	return out
end

function f()
	L = 6
	N = L
	state = onesState(L)
	#measOp = singleSubspaceProjectors(L, N)
	measOp = generateProjectionOperators(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=10, measOp=measOp, traj=100, dt=0.02, time=10.0, p=0.0, f=1.0, U=0.14, J=1.0, Ψ₀=state)
	probabilities = [0.01, 0.06, 0.08, 0.15]
	results = calcWithDifferentProb(p, probabilities)
	pl = plot(p.t.times, results[1], label="p = $(probabilities[1])", ylim=(0.0, 2.5))
	for i in 2:length(results)
		plot!(p.t.times, results[i], label="p = $(probabilities[i])")
	end
	return pl
end

f()

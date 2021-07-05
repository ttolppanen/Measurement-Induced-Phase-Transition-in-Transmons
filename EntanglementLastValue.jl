using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function calcWithDifferentProb(p::Parameters, probabilities)
	out = []
	outVar = []
	for prob in probabilities
		param = ParametersConstructorWithP(p, prob)
		@time sol = MIPT(param, onlyLastValue = true, projectAfterTimeStep = true)
		@time res, var = calcMeanAndVar(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, Int(p.L / 2)))
		push!(out, res[1])
		push!(outVar, var[1])
	end
	return out, outVar
end
function entanglementAsAFunctionOfP(L, probabilities)
	N = L
	state = onesState(L)
	#measOp = singleSubspaceProjectors(L, N)
	measOp = generateProjectionOperators(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=10, measOp=measOp, traj=10, dt=0.02, time=10.0, p=0.0, f=1.0, U=0.14, J=1.0, Ψ₀=state)
	return calcWithDifferentProb(p, probabilities)
end

function f()
	probabilities = 0.02:0.02:0.1
	L = 8
	results, var = entanglementAsAFunctionOfP(L, probabilities)
	plot!(probabilities, results, ribbon=var, label="L = $L", ylim=(0.0, 2.5))
end

f()

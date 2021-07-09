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
function entanglementAsAFunctionOfP(L, probabilities, traj)
	N = floor(Int, L/2)
	state = zeroOneState(L, N)
	#measOp = singleSubspaceProjectors(L, N)
	measOp = generateProjectionOperators(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=10, measOp=measOp, traj=traj, dt=0.02, time=10.0, p=0.0, f=1.0, U=4.0, J=1.0, Ψ₀=state)
	return calcWithDifferentProb(p, probabilities)
end

function f(L, traj)
	probabilities = 0.02:0.02:0.14
	results, var = entanglementAsAFunctionOfP(L, probabilities, traj)
	plot(probabilities, results, ribbon=sqrt.(var), fillalpha=0.15, label="L = $L")
end
function f!(L, traj)
	probabilities = 0.02:0.02:0.14
	results, var = entanglementAsAFunctionOfP(L, probabilities, traj)
	plot!(probabilities, results, ribbon=sqrt.(var), fillalpha=0.15, label="L = $L")
end

f(4, 10000)
f!(6, 6000)
f!(8, 3000)

using Plots, MIPTM, Distributions
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function calcWithDifferentProb(p::Parameters, probabilities)
	out = []
	outVar = []
	for prob in probabilities
		param = ParametersConstructorWithP(p, prob)
		@time sol = MIPT(param, onlyLastValue = true, projectAfterTimeStep = true)
		#@time res, var = calcMeanAndVar(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, Int(p.L / 2), p.cap))
		@time res, var = calcMeanAndVar(sol, Ψ->halfBosonNumber(Ψ, p.L, p.N, p.cap))
		push!(out, res[1])
		push!(outVar, var[1])
	end
	return out, outVar
end
function makeParam(L, traj)
	#N = floor(Int, L/2)
	N = L
	state = onesState(L, N)
	measOp = singleSubspaceProjectors(L, N)
	#measOp = generateProjectionOperators(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=3, measOp=measOp, traj=traj, dt=0.02,
	time=30.0, p=0.0, f=1.0, U=0.14, J=1.0, Ψ₀=state, mean=0.0, stantardDeviation=1.0)

	return p
end

function f(L, traj)
	probabilities = 0.02:0.01:0.09
	p = makeParam(L[1], traj[1])
	results = []
	var = []
	res1, var1 = calcWithDifferentProb(p, probabilities)
	push!(results, res1)
	push!(var, var1)
	pl = plot(probabilities, results[1], ribbon=sqrt.(var[1]), fillalpha=0.15, label="L = $(L[1])")
	for i in 2:length(traj)
		p = makeParam(L[i], traj[i])
		resᵢ, varᵢ = calcWithDifferentProb(p, probabilities)
		push!(results, resᵢ)
		push!(var, varᵢ)
		plot!(probabilities, results[i], ribbon=sqrt.(var[i]), fillalpha=0.15, label="L = $(L[i])")
	end
	savePlotData(probabilities, (results, var), "Fluc_S_2000_1000_100", p, "1111..."; notes="")
	display(pl)
end

f([4, 6, 8], [2000, 1000, 100])

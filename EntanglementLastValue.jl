using Plots, MIPTM, Distributions
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function calcWithDifferentProb(p::Parameters, probabilities)
	out1 = []
	outVar1 = []
	out2 = []
	outVar2 = []
	for prob in probabilities
		param = ParametersConstructorWithP(p, prob)
		@time sol = MIPT(param, onlyLastValue = true, projectAfterTimeStep = true)
		@time res1, var1 = calcMeanAndVar(sol, Ψ->entanglement_entropy(p.L, p.N, Ψ, Int(p.L / 2), p.cap))
		@time res2, var2 = calcMeanAndVar(sol, Ψ->halfBosonNumber(Ψ, p.L, p.N, p.cap))
		push!(out1, res1[1])
		push!(outVar1, var1[1])
		push!(out2, res2[1])
		push!(outVar2, var2[1])
	end
	return out1, outVar1, out2, outVar2
end
function makeParam(L, traj)
	#N = floor(Int, L/2)
	N = Int(L/2)
	cap = 1
	state = zeroOneState(L, N, cap)
	#measOp = singleSubspaceProjectors(L, N, cap)
	measOp = generateProjectionOperators(L, N, cap)
	p = ParametersConstructor(L=L, N=N, cap=cap, sdim=3, measOp=measOp, traj=traj, dt=0.02,
	time=30.0, p=0.0, f=1.0, U=0.14, J=1.0, Ψ₀=state, mean=0.0, stantardDeviation=0.0)

	return p
end
function makePlot(probabilities, res, var, L)
	pl = plot(probabilities, res[1], ribbon=sqrt.(var[1]), fillalpha=0.15, label="L = $(L[1])")
	for i in 2:length(res)
		plot!(probabilities, res[i], ribbon=sqrt.(var[i]), fillalpha=0.15, label="L = $(L[i])")
	end
	return pl
end

function f(L, traj)
	probabilities = 0.02:0.01:0.09
	p = makeParam(L[1], traj[1])
	results1 = []
	variances1 = []
	results2 = []
	variances2 = []
	res1, var1, res2, var2 = calcWithDifferentProb(p, probabilities)
	push!(results1, res1)
	push!(variances1, var1)
	push!(results2, res2)
	push!(variances2, var2)
	for i in 2:length(traj)
		p = makeParam(L[i], traj[i])
		res1ᵢ, var1ᵢ, res2ᵢ, var2ᵢ = calcWithDifferentProb(p, probabilities)
		push!(results1, res1ᵢ)
		push!(variances1, var1ᵢ)
		push!(results2, res2ᵢ)
		push!(variances2, var2ᵢ)
	end
	makePlot(probabilities, results1, variances1, L)
	savePlotData(probabilities, (results1, variances1), "ELV_S_Cap1_State0101_10000_10000_10000", p, "0101..."; notes="Total space projection")
	makePlot(probabilities, results2, variances2, L)
	savePlotData(probabilities, (results2, variances2), "Fluc_S_Cap1_State0101_10000_10000_10000", p, "0101..."; notes="Total space projections")
end

f([4, 6, 8], [10000, 10000, 10000])

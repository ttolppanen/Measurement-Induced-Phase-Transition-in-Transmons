using Plots, MIPTM, Distributions, ParametersModule
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function calcWithDifferentProb(p::Parameters, probabilities)
	out = []
	outVar = []
	ent(Ψ) = entanglement_entropy(p.sp.L, p.sp.N, Ψ, Int(p.sp.L / 2), p.cap)
	fluc(Ψ) = halfBosonNumber(Ψ, p.sp.L, p.sp.N, p.cap)
	functions = [ent, fluc]
	for _ in 1:length(functions)
		push!(out, [])
		push!(outVar, [])
	end
	push!(out, [])
	push!(outVar, [])

	for prob in probabilities
		p.pp.p = prob
		@time sol = MIPT(p, onlyLastValue = true, projectAfterTimeStep = true)
		for i in 1:length(functions)
			@time res, var = calcMeanAndVar(sol, functions[i])
			push!(out[i], res[1])
			push!(outVar[i], var[1])
		end
		@time res = properFluc(sol, p)
		push!(out[3], res[1])
	end
	return out, outVar
end
function makeParam(L, traj)
	#N = floor(Int, L/2)
	N = L
	cap = 2
	state = onesState(L, cap)
	projOp = singleSubspaceProjectors(L, N, cap)
	#projOp = generateProjectionOperators(L, N, cap)
	sp = SystemParameters(L=L, N=N, Ψ₀=state)
	pp = ProjectionParameters(p=0.0, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters(L=L, N=N, cap=cap, U=-0.14, Uσ=0.07, w=0.0, wσ = 0.015)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, cap=cap, sdim=3, dt=0.02, time=30.0, traj=traj)
	return p
end
function makePlot(probabilities, res, var, L)
	pl = plot(probabilities, res[1], ribbon=sqrt.(var[1]), fillalpha=0.15, label="L = $(L[1])")
	for i in 2:length(res)
		plot!(probabilities, res[i], ribbon=sqrt.(var[i]), fillalpha=0.15, label="L = $(L[i])")
	end
	return pl
end
function makePlot(probabilities, res, L)
	pl = plot(probabilities, res[1], fillalpha=0.15, label="L = $(L[1])")
	for i in 2:length(res)
		plot!(probabilities, res[i], fillalpha=0.15, label="L = $(L[i])")
	end
	return pl
end

function f(L, traj)
	probabilities = 0.0:0.01:0.09
	results = []
	variances = []
	numOfFunctions = 3
	for _ in 1:numOfFunctions
		push!(results, [])
		push!(variances, [])
	end
	p = makeParam(L[1], traj[1])
	for i in 1:length(L)
		p = makeParam(L[i], traj[i])
		res, var = calcWithDifferentProb(p, probabilities)
		for j in 1:numOfFunctions
			push!(results[j], res[j])
			push!(variances[j], var[j])
		end
	end
	pl = makePlot(probabilities, results[1], variances[1], L)
	savePlotData(probabilities, (results[1], variances[1]), "ELV_S_d3_dis_Attractive_Superfluid_7000_2000_100", p, "1111..."; notes="|1><1| projection")
	makePlot(probabilities, results[2], variances[2], L)
	savePlotData(probabilities, (results[2], variances[2]), "Fluc_S_d3_dis_Attractive_Superfluid_7000_2000_100", p, "1111..."; notes="|1><1| projection")
	makePlot(probabilities, results[3], L)
	savePlotData(probabilities, (results[3], results[3]), "ProperFluc_S_d3_dis_Attractive_Superfluid_7000_2000_100", p, "1111..."; notes="|1><1| projection, variances are not real data! They are just the same as the result")
	display(pl)
end

f([4, 6, 8], [7000, 2000, 100])

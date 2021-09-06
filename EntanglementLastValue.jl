using Plots, MIPTM, Distributions, ParametersModule
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function calcWithDifferentProb(p::Parameters, probabilities)
	out = []
	outVar = []
	ent(Ψ) = entanglement_entropy(p.sp.L, p.sp.N, Ψ, Int(p.sp.L / 2), cap=p.sp.cap)
	fluc(Ψ) = halfBosonNumber(Ψ, p.sp.L, p.sp.N, cap=p.sp.cap)
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
		display("Result times")
		for i in 1:length(functions)
			@time res, var = calcMeanAndVar(sol, functions[i])
			push!(out[i], res[1])
			push!(outVar[i], var[1])
		end
		@time res = properFluc(sol, p)
		display("end results")
		push!(out[end], res[1])
	end
	return out, outVar
end
function makeParam(L, traj)
	#N = floor(Int, L/2)
	N = L
	sp = SystemParameters(L=L, N=N, cap=2)
	state = onesState(sp)
	projOp = singleSubspaceProjectors(sp)
	#projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=0.0, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters(sp=sp, U=-5.14, Uσ=0.07, w=0.0, wσ=0.015)
	#bhp = BoseHubbardParameters(sp=sp, U=-5.14)
	#bhp = BoseHubbardParameters(sp=sp, U=0.0, J=1.0, Jσ=0.5)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=3, dt=0.02, time=30.0, traj=traj)
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
	#plot!(probabilities, x->2.0/((x/0.02)^2 + 3))
	savePlotData(probabilities, (results[1], variances[1]), "Fixed_ELV_S_d3_dis_Attractive_Insulator_10000_5000_300", p, "1111..."; notes="|1><1| projection")
	makePlot(probabilities, results[2], variances[2], L)
	savePlotData(probabilities, (results[2], variances[2]), "Fixed_Fluc_S_d3_dis_Attractive_Insulator_10000_5000_300", p, "1111..."; notes="|1><1| projection")
	makePlot(probabilities, results[3], L)
	savePlotData(probabilities, (results[3], results[3]), "Fixed_ProperFluc_S_d3_dis_Attractive_Insulator_10000_5000_300", p, "1111..."; notes="|1><1| projection, variances are not real data! They are just the same as the result")
	display(pl)
	#pl = makePlot(probabilities, results[1], L)
	#return pl
end

f([4, 6, 8], [10000, 5000, 300])
#f([2], [5000])

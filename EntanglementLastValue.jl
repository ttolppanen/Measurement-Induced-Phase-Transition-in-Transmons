using Plots, MIPTM, Distributions, ParametersModule
using BSON: @save, @load
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function convertResults(out, outVar, functions)

end
function calcWithDifferentProb(p::Parameters, probabilities)
	out = []
	outVar = []
	ent(Ψ) = entanglement_entropy(p.sp.L, p.sp.N, Ψ, Int(p.sp.L / 2), cap=p.sp.cap)
	fluc(Ψ) = halfBosonNumber(Ψ, p.sp.L, p.sp.N, p.sp.dim, cap=p.sp.cap)
	functions = [fluc]
	for _ in 1:length(functions)
		push!(out, [])
		push!(outVar, [])
	end
	#push!(out, [])
	#push!(outVar, [])
	for prob in probabilities
		p.pp.p = prob
		sol = MIPT(p, onlyLastValue = true, projectAfterTimeStep = true)
		for i in 1:length(functions)
			res, var = calcMeanAndVar(sol, functions[i])
			push!(out[i], res[1])
			push!(outVar[i], var[1])
		end
		#res = properFluc(sol, p)
		#display("end results")
		#push!(out[end], res[1])
	end
	return out, outVar
end
function makeParam(L, traj)
	#N = floor(Int, L/2)
	N = L
	sp = SystemParameters(L=L, N=N)
	state = onesState(sp)
	projOp = singleSubspaceProjectors(sp)
	#projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=0.0, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters(sp=sp, U=0.14)
	#bhp = BoseHubbardParameters(sp=sp, U=-5.14)
	#bhp = BoseHubbardParameters(sp=sp, U=0.0, J=1.0, Jσ=0.5)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=3, dt=0.02, time=30.0, traj=traj)
	return p
end
function solveWithDifferentSizes(sizes::Array{Tuple{Int64, Int64}, 1}, functions)
	results = Dict{Int64, Result}()
	L = sizes[1][1]
	traj = sizes[1][2]
	p = makeParam(L, traj)
	for i in 1:length(L)
		if i != 1
			L = sizes[i][1]
			traj = sizes[i][2]
			p = makeParam(L, traj)
		end
		res, var = calcWithDifferentProb(p, probabilities)
		for j in 1:Length(functions)
			push!(results[j], res[j])
			push!(variances[j], var[j])
		end
	end
	return results, p
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
function replaceData(title, index, newData)
	path = pwd() * "/Plots/" * title * "/data.bson"
	@load path x y
	y[1][index] = newData[1]
	y[2][index] = newData[2]
	@save path x y
end
function addToData(title, newData)
	path = pwd() * "/Plots/" * title * "/data.bson"
	@load path x y
	push!(y[1], newData[1])
	push!(y[2], newData[2])
	@save path x y
end
struct Result
	mean
	variance
	function Result(res::Tuple{Any, Any})
		new(res[1], res[2])
	end
	function Result(res::Array{<:Any, 1})
		new(res, nothing)
	end
end

function f(L, traj)
	probabilities = [0.01, 0.02, 0.03, 0.04, 0.05, 0.055, 0.06, 0.0625, 0.065, 0.0675, 0.07, 0.075, 0.08, 0.09]
	#probabilities = collect(0.01:0.02:0.09)
	#probabilities .*= 5
	results = []
	variances = []
	numOfFunctions = 1
	for _ in 1:numOfFunctions
		push!(results, [])
		push!(variances, [])
	end
	p = makeParam(L[1], traj[1])
	for i in 1:length(L)
		p = makeParam(L[i], traj[i])
		display("res and var time")
		@time res, var = calcWithDifferentProb(p, probabilities)
		for j in 1:numOfFunctions
			push!(results[j], res[j])
			push!(variances[j], var[j])
		end
	end
	pl = makePlot(probabilities, results[1], variances[1], L)
	#replaceData("ELV_S_5000_1000_100", 3, (results[1][1], variances[1][1]))
	#addToData("ELV_S_5000_1000_300", (results[1], variances[1]))
	#plot!(probabilities, x->2.0/((x/0.02)^2 + 3))
	savePlotData(probabilities, (results[1], variances[1]), "Fluc_S_5000_1000_300", p, "1111..."; notes="|1><1| projection")
	#makePlot(probabilities, results[2], variances[2], L)
	#savePlotData(probabilities, (results[2], variances[2]), "202020test_5000_2000_100", p, "2020..."; notes="|1><1| projection")
	#makePlot(probabilities, results[3], L)
	#savePlotData(probabilities, (results[3], results[3]), "202020test_fluc_5000_2000_100", p, "2020..."; notes="|1><1| projection")
	display(pl)
	#pl = makePlot(probabilities, results[1], L)
	#return pl
end

f([4, 6, 8], [5000, 1000, 300])
#f([2], [5000])

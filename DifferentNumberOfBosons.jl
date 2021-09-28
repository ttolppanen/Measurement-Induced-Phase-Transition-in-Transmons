using MIPTM, ParametersModule, Plots
include.(["OllisCode/Entropy.jl", "OllisCode/Basis.jl"])

function calcWithDifferentProb(probabilities, p::Parameters, funcs)
	out = []
	for _ in 1:length(funcs)
		push!(out, Result())
	end
	for prob in probabilities
		p.pp.p = prob
		sol = MIPT(p, onlyLastValue = true, projectAfterTimeStep = true)
		for i in 1:length(funcs)
			newRes = funcs[i](sol, p)
			addResult!(out[i], newRes)
		end
	end
	return out
end
function calcAll(prob, sizesAndTraj, funcs)
	out = []
	for st in sizesAndTraj
		p = makeParam(st)
		res = calcWithDifferentProb(prob, p, funcs)
		push!(out, res)
	end
	return out
end
struct Result
	mean
	variance
	function Result()
		new([], [])
	end
end
function addResult!(result::Result, newRes::Tuple{Float64, Float64})
	push!(result.mean, newRes[1])
	push!(result.variance, newRes[2])
end
function addResult!(result::Result, newRes::Float64)
	push!(result.mean, newRes)
end

struct ST #SizeAndTrajectory
	size::Int64
	traj::Int64
end
function convertResults(sol)
	numOfSizes = length(sol)
	numOfFunctions = length(sol[1])
	out = []
	for i in 1:numOfFunctions
		push!(out, [])
		for j in 1:numOfSizes
			push!(out[end], sol[j][i])
		end
	end
	return out
end
function convertForSize(sol)
	out = []
	for i in sol
		push!(out, i[1])
	end
	return out
end
function plotForFunction(probabilities, results, i)
	numOfSizes = length(results[1])
	if length(results[i][1].variance) == 0
		pl = plot(probabilities, results[i][1].mean)
		for j in 2:numOfSizes
			plot!(probabilities, results[i][j].mean)
		end
		return pl
	else
		pl = plot(probabilities, results[i][1].mean, ribbon=sqrt.(results[i][1].variance), fillalpha=0.15)
		for j in 2:numOfSizes
			plot!(probabilities, results[i][j].mean, ribbon=sqrt.(results[i][j].variance), fillalpha=0.15)
		end
		return pl
	end
end
function functionsToCalculate()
	entForState(Ψ, p) = entanglement_entropy(p.sp.L, p.sp.N, Ψ, Int(p.sp.L / 2), cap=p.sp.cap)
	entanglement(sol, p) = calcMeanAndVar(sol, Ψ -> entForState(Ψ, p))
	return [entanglement]
end
function makeParam(st::ST)
	sp = SystemParameters(L=6, N=st.size)
	state = firstState(sp)
	#projOp = singleSubspaceProjectors(sp)
	projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=0.0, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters(sp=sp, U=5.0)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=3, dt=0.02, time=30.0, traj=st.traj)
end
function firstState(sp)
	L = sp.L; N = sp.N; dim = sp.dim; cap = sp.cap
	basisState::Array{Int64,1} = first_state(L, N, cap)
	state = zeros(dim)
	state[find_index(basisState, cap)] = 1.
	return state
end

function f()
	probabilities = [0.08]
	sizesAndTraj = []
	maxN = 8
	for i in 2:maxN
		push!(sizesAndTraj, ST(i, 500))
	end
	funcs = functionsToCalculate()
	results = calcAll(probabilities, sizesAndTraj, funcs) #size[functions[result1, result2,...],...]

	#=
	p = makeParam(sizesAndTraj[1])
	results = convertResults(results)
	pl = plotForFunction(probabilities, results, 1)
	mean = [i.mean for i in results[1]]
	var = [i.variance for i in results[1]]
	#savePlotData(probabilities, (mean, var), "NTestU5", p, "|N000000>"; notes="Total Projection")
	=#

	results = convertForSize(results)
	sizes = [st.size for st in sizesAndTraj]
	mean = [i.mean[1] for i in results]
	var = [i.variance[1] for i in results]
	pl = plot(sizes, mean, ribbon=sqrt.(var))
	display(pl)
end

f()

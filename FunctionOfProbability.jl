using MIPTM, ParametersModule, Plots
include.("OllisCode/Entropy.jl")

function calcWithDifferentProb(probabilities, p::Parameters, funcs)
	out = []
	for _ in 1:length(funcs)
		push!(out, Result())
	end
	for prob in probabilities
		p.pp.p = prob
		@time sol = MIPT(p, onlyLastValue = true, projectAfterTimeStep = true)
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
	L::Int64
	N::Int64
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
	state(sol, p) = calcMeanAndVar(sol, Ψ -> expVal(Ψ, firstState(p.sp) * firstState(p.sp)'))
	return [entanglement] #These need to be functions that take in the solution (sol) and the parameters (p), so f(sol, p)
end
function makeParam(st::ST)
	sp = SystemParameters(L=st.L, N=st.N)
	state = onesState(sp)
	projOp = singleSubspaceProjectors(sp)
	#projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=0.0, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters(sp=sp, U=5.0)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=3, dt=0.02, time=30.0, traj=st.traj)
end
function meanAndVar(results)
	return [res.mean for res in results], [res.variance for res in results]
end
function firstState(sp)
	L = sp.L; N = sp.N; dim = sp.dim; cap = sp.cap
	basisState::Array{Int64,1} = first_state(L, N, cap)
	state = zeros(dim)
	state[find_index(basisState, cap)] = 1.
	return state
end

function f()
	probabilities = 0.0:0.001:0.01
	sizesAndTraj = [ST(4,4, 1000), ST(6,6, 500), ST(8,8, 60)]
	funcs = functionsToCalculate()
	results = calcAll(probabilities, sizesAndTraj, funcs) #size[functions[result1, result2,...],...]

	p = makeParam(sizesAndTraj[1])
	results = convertResults(results) #functions[size[result1, result2,...],...]
	pl = plotForFunction(probabilities, results, 1)
	#savePlotData(probabilities, meanAndVar(results[1]), "DisorderTest", p, "|N00000>"; notes="Total Projection")
	#plotForFunction(probabilities, results, 2)
	#savePlotData(probabilities, meanAndVar(results[2]), "Fluc_N_L_15", p, "|N00000>"; notes="Total Projection")
	display(pl)
end

f()

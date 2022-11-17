using MIPTM, ParametersModule, Plots
include.("OllisCode/Entropy.jl")

function calcWithDifferentProb(probabilities, p::Parameters, funcs)
	out = []
	for _ in 1:length(funcs)
		push!(out, [])
	end
	for prob in probabilities
		p.pp.p = prob
		@time sol = MIPT(p, onlyLastValue = true, projectAfterTimeStep = true)
		for i in eachindex(funcs)
			newRes = funcs[i](sol, p)
			push!(out[i], newRes)
		end
	end
	return out
end
function calcAll(prob, sizesAndTraj, funcs)
	out = []
	for st in sizesAndTraj
		p = makeParam(st)
		println("L = $(p.sp.L)")
		@time res = calcWithDifferentProb(prob, p, funcs)
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

function halfL_BosonNumbermean(Ψ, L) # Numbers of bosons in the half of the chain
	nₕ = number(L, 1)
	for i in 2:Int(round(L/2))
		nₕ .+= number(L, i)
	end
	return expVal4(Ψ, nₕ)
end
function L_BosonNumbermean(Ψ, L) # Number of boson in the total chain
	nₕ = number(L, 1)
	for i in 2:L
		nₕ .+= number(L, i)
	end
	return expVal4(Ψ, nₕ)
end
function l_BosonNumbermean(Ψ, L) # Number of bosons at l in the middle of the chain
	nₕ = number(L, round(Int64,L/2))
	#nₕ = number(L, L)
	return expVal4(Ψ, nₕ)
end

function corrl_BosonNumbermean(Ψ, L) # Correlation of the number of bosons at l and l+1
	nₕa = number(L, round(Int64,L/2))
	nₕb = number(L, round(Int64,L/2+1))
	return expVal4(Ψ, nₕa.*nₕb)
end

#2) Observables to be computed as real measurement with the function real_measurement_results
function halfL_BosonNumber(L) # OPERATOR for the number of bosons in half of the chain
	nₕ = number(L, 1)
	for i in 2:round(Int64,L/2)
		nₕ .+= number(L, i)
	end
	return nₕ
end
function L_BosonNumber(L) # OPERATOR for the total number of bosons
	nₕ = number(L, 1)
	for i in 2:L
		nₕ .+= number(L, i)
	end
	return nₕ
end
function corrl_BosonNumber(L) # OPERATOR for the correlation of the number of bosons at l and l+1
	nₕa = number(L, round(Int64,L/2))
	nₕb = number(L, round(Int64,L/2+1))
	return nₕa.*nₕb
end

function functionsToCalculate()
	#these are the expectation values (or postselected) observables mean values
	entForState(Ψ, p) = entanglement_entropy(p.sp.L, Ψ, Int(p.sp.L / 2))
	entanglement(sol, p) = calcMeanAndVar4(sol, Ψ -> entForState(Ψ, p))
	fluctuations(sol, p) = calcMeanAndVar4(sol, Ψ -> Fluc4(Ψ, p.sp.L))

	halfLNexpval(sol, p) = calcMeanAndVar4(sol, Ψ -> halfL_BosonNumbermean(Ψ, p.sp.L))
	LNexpval(sol, p) = calcMeanAndVar4(sol, Ψ -> L_BosonNumbermean(Ψ, p.sp.L))
	lNexpval(sol, p) = calcMeanAndVar4(sol, Ψ -> l_BosonNumbermean(Ψ, p.sp.L))
	corrlNexpval(sol, p) = calcMeanAndVar4(sol, Ψ -> corrl_BosonNumbermean(Ψ, p.sp.L))

	#WE NEED TO INCLUDE HERE SOMETHING LIKE THIS
	halfLN(sol,p)=real_measurement_results(sol, halfL_BosonNumber(p.sp.L))
	LN(sol,p)=real_measurement_results(sol, L_BosonNumber(p.sp.L))
	lN(sol, p)=real_measurement_results(sol, number(p.sp.L, round(Int64,p.sp.L/2)))
	corrlN(sol,p)=real_measurement_results(sol, corrl_BosonNumber(p.sp.L))

	#return [entanglement,properFluc4,fluctuations,halfLNexpval,LNexpval,lNexpval,corrlNexpval]
	return [entanglement,properFluc4,fluctuations,halfLNexpval,LNexpval,lNexpval,corrlNexpval,halfLN,LN,lN,corrlN] #THIS IS THE PROPER ONE
end
function makeParam(st::ST)
	sp = SystemParameters(L=st.L, N=st.N)
	state = onesZeroState(sp)
	#projOp = singleSubspaceProjectors(sp)
	projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=0.0, f=1.0, projOp=projOp)
	bhp = BoseHubbardParameters4(sp=sp, w=100.0, wσ=0.0,U=10.0, Uσ=0.0,J=1.0, Jσ=0.0) #Minus sign for U is included in HU, U>0 here is attractive
	return Parameters4(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=6, dt=0.02, time=20.0, traj=st.traj)#sdim is Krylov subspace dim
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

function f(measurementType)
	#SELECT PROBABILITY RANGE AND TRAJECTORIES
	probabilities = 0.005:0.001:0.1
	#probabilities = 0.1:0.025:1.0
	traj = 10
	sizesAndTraj = [ST((4,2), traj), ST((6,3), traj), ST((8,4), traj), ST((10,5), traj)]
	#sizesAndTraj = [ST((12,6), traj)]

	funcs = functionsToCalculate()
	results = calcAll(probabilities, sizesAndTraj, funcs, measurementType) #size[functions[result1, result2,...],...]

	p = makeParam(sizesAndTraj[1], measurementType)
	results = convertResults(results) #functions[size[result1, result2,...],...]

	# HERE WE SAVE ALL THE CALCULATED FUNCTIONS
	folder = "Data"
	notes = "Measurement_" * string(measurementType)
	savetokataja = true
	savePlotData4(probabilities, results[1], "M$(measurementType)_Entanglement", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #THE NOTES COULD INCLUDE THE TYPE OF THE MEASUREMENT (1,2,3)
	savePlotData4(probabilities, results[2], "M$(measurementType)_Proper_fluctuation", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, results[3], "M$(measurementType)_fluctuations", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, results[4], "M$(measurementType)_halfLNexpval", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, results[5], "M$(measurementType)_LNexpval", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, results[6], "M$(measurementType)_smallLNexpval", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, results[7], "M$(measurementType)_corrlNexpval", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	#data where the superposition is destroyed
	savePlotData4(probabilities, results[8], "M$(measurementType)_halfLN", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, results[9], "M$(measurementType)_LN", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, results[10], "M$(measurementType)_smallLN", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, results[11], "M$(measurementType)_corrlN", p, "|1010..>"; notes=notes, folder=folder, savetokataja=savetokataja) #results[1] ent, results[2] ProperFluc....
	
	#=
	# THESE NEEDS TO BE INCLUDED IN THE SAVED FILES
	savePlotData4(probabilities, meanAndVar(results[8]), "halfLN", p, "|1010..>"; notes="Measurement_1") #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, meanAndVar(results[9]), "LN", p, "|1010..>"; notes="Measurement_1") #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, meanAndVar(results[10]), "lN", p, "|1010..>"; notes="Measurement_1") #results[1] ent, results[2] ProperFluc....
	savePlotData4(probabilities, meanAndVar(results[11]), "corrlN", p, "|1010..>"; notes="Measurement_1") #results[1] ent, results[2] ProperFluc....

	#THESE PLOTS ARE JUST FOR CHECKING PURPOSES SINCE THE FINAL PLOTS ARE GOING TO BE DONE ELSEWHERE
	pl1 = plotForFunction(probabilities, results, 1, sizesAndTraj)
	pl2 = plotForFunction(probabilities, results, 2, sizesAndTraj)
	pl3 = plotForFunction(probabilities, results, 3, sizesAndTraj)
	pl4 = plotForFunction(probabilities, results, 4, sizesAndTraj)
	pl5 = plotForFunction(probabilities, results, 5, sizesAndTraj)
	pl6 = plotForFunction(probabilities, results, 6, sizesAndTraj)
	pl7 = plotForFunction(probabilities, results, 7, sizesAndTraj)
	#pl8 = plotForFunction(probabilities, results, 8, sizesAndTraj)
	#pl9 = plotForFunction(probabilities, results, 9, sizesAndTraj)
	#pl10 = plotForFunction(probabilities, results, 10, sizesAndTraj)
	#pl11 = plotForFunction(probabilities, results, 11, sizesAndTraj)

	display(pl1)
	display(pl2)
	display(pl3)
	display(pl4)
	display(pl5)
	display(pl6)
	display(pl7)
	#display(pl8)
	#display(pl9)
	#display(pl10)
	#display(pl11)
	=#

end

f()

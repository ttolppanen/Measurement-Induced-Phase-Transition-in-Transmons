module MIPTM
	using DifferentialEquations, IterTools, LinearAlgebra, SparseArrays, Plots
	using Distributions, ParametersModule
	using Statistics: mean
	using BSON: @save
	include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl", "OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

	export calcMean, calcMeanAndVar, expVal
	export MIPT
	export singleSubspaceProjectors, onesState, zeroOneState, oneZeroState, twoZeroState
	export generateProjectionOperators, generateSingleSite
	export halfBosonNumber, properFluc
	export savePlotData, readBsonFile

	function generateProjectionOperators(L, N, dim; cap=N)
		out = []
		for l in 1:L
			oneSiteOperators = []
			for n in 0:min(N, cap)
				push!(oneSiteOperators, projector(L, N, dim, l, n, cap=cap))
			end
			push!(out, oneSiteOperators)
		end
		return out
	end
	function generateProjectionOperators(sp::SystemParameters)
		return generateProjectionOperators(sp.L, sp.N, sp.dim, cap=sp.cap)
	end
	function generateSingleSite(L, N, f::Function)
		out = []
		for l in 1:L
			oneSiteOperators = []
			push!(oneSiteOperators, f(L, N, l))
			push!(out, oneSiteOperators)
		end
		return out
	end
	function singleSubspaceProjectors(L, N, dim; cap=N)
		f(L, N, l) = projector(L, N, dim, l, 1, cap=cap)
		return generateSingleSite(L, N, f)
	end
	function singleSubspaceProjectors(sp::SystemParameters)
		singleSubspaceProjectors(sp.L, sp.N, sp.dim, cap=sp.cap)
	end
	function halfBosonNumber(Ψ, L, N, dim; cap=N) #If this is slow change number to take dim straigth...
		nₕ = number(L, N, dim, 1, cap=cap)
		for i in 2:Int(round(L/2))
			nₕ .+= number(L, N, dim, i, cap=cap)
		end
		return expVal(Ψ, nₕ^2) - expVal(Ψ, nₕ)^2
	end
	function onesState(L::Int64, dim; cap=L)
		state = zeros(dim)
		state[find_index(ones(Int64, L), cap)] = 1.
		return state
	end
	function onesState(sp::SystemParameters)
		return onesState(sp.L, sp.dim, cap=sp.cap)
	end
	function oneZeroState(L::Int64, N::Int64, dim; cap=N)
		basisState::Array{Int64,1} = []
		for i in 1:L
			push!(basisState, i%2)
		end
		state = zeros(dim)
		state[find_index(basisState, cap)] = 1.
		return state
	end
	function twoZeroState(L::Int64, N::Int64, dim; cap=N)
		basisState::Array{Int64,1} = []
		for i in 1:L
			push!(basisState, 2*(i%2))
		end
		state = zeros(dim)
		state[find_index(basisState, cap)] = 1.
		return state
	end
	function twoZeroState(sp::SystemParameters)
		return twoZeroState(sp.L, sp.N, sp.dim, cap=sp.cap)
	end
	function oneZeroState(sp::SystemParameters)
		return oneZeroState(sp.L, sp.N, sp.dim, cap=sp.cap)
	end
	function zeroOneState(L::Int64, N::Int64, dim; cap=N)
		basisState::Array{Int64,1} = [0]
		for i in 1:L-1
			push!(basisState, i%2)
		end
		state = zeros(dim)
		state[find_index(basisState, cap)] = 1.
		return state
	end
	function zeroOneState(sp::SystemParameters)
		return zeroOneState(sp.L, sp.N, sp.dim, cap=sp.cap)
	end
	function expVal(s::Array{Complex{Float64},1}, op::Union{Array{Float64,2}, Array{Complex{Float64},2}, SparseMatrixCSC{Float64,Int64}})#Jos s on ket
		return real(s' * op * s)
	end
	function expVal(ρ::Array{Complex{Float64},2}, op::Array{Complex{Float64},2})#Tiheysoperaattorille
		return real(tr(op*ρ))
	end
	function expVal(ρ::Array{Complex{Float64},2}, op::Array{Complex{Float64},2}, mPA1::Array{Complex{Float64},2})#Tiheysoperaattorille
		mul!(mPA1, op, ρ)
		return real(tr(mPA1))
	end
	function calcMean(sol, f::Function)
		numOfVal = length(sol)

		threads = []
		for _ in 1:Threads.nthreads()
			push!(threads, zeros(length(sol[1])))
		end

        Threads.@threads for i in 1:numOfVal
            threads[Threads.threadid()] .+= f.(sol[i])
        end
		mean = sum(threads)
		mean ./= numOfVal
		if length(mean) == 1
			return mean[1]
		end
		return mean
    end
	function calcMeanAndVar(sol, f::Function)
		numOfVal = length(sol)

		threadsMean = []
		threadsVar = []
		for _ in 1:Threads.nthreads()
			push!(threadsMean, zeros(length(sol[1])))
			push!(threadsVar, zeros(length(sol[1])))
		end

		Threads.@threads for i in 1:numOfVal
			fVal = f.(sol[i])
            threadsMean[Threads.threadid()] .+= fVal
			threadsVar[Threads.threadid()] .+= fVal.^2
        end
		mean = sum(threadsMean)
		var = sum(threadsVar)
		mean .= mean./numOfVal
		var .= var./numOfVal .- mean.^2
		var .= round.(var, digits=12)
		if length(mean) == 1
			return mean[1], var[1]
		end
        return mean, var
    end
	function properFluc(sol, p::Parameters)
		nₕ = number(p.sp.L, p.sp.N, p.sp.dim, 1, cap=p.sp.cap)
		for i in 2:Int(round(p.sp.L/2))
			nₕ .+= number(p.sp.L, p.sp.N, p.sp.dim, i, cap=p.sp.cap)
		end
		f1(Ψ) = expVal(Ψ, nₕ)
		f2(Ψ) = expVal(Ψ, nₕ.^2)
		d = calcMean(sol, f1)
		return calcMean(sol, f2) - calcMean(sol, f1).^2
    end
	function measurementEffect!(Ψ, p::Parameters)
		for l in 1:p.sp.L
			if rand(Float64) < p.pp.p  #Does the measurement happen?
				probForProjection = rand(Float64)
				pⱼ = 0 #Probability for a single projection
				for j in 1:length(p.pp.projOp[l])
					if j == length(p.pp.projOp[l])
						projection!(Ψ, p.pp.projOp[l][j])
						break
					else
						pⱼ += projectionProbability(Ψ, p.pp.projOp[l][j])
						if probForProjection < pⱼ
							projection!(Ψ, p.pp.projOp[l][j])
							break
						end
					end
				end
			end
		end
	end
	function projection!(Ψ, op)
		Ψ .= op * Ψ
		Ψ ./= norm(Ψ)
	end
	function projectionProbability(Ψ, op)
		return expVal(Ψ, op' * op)
	end
	function lastValues(sol)
		res = []
		for i in sol
			push!(res, [last(i)])
		end
		return res
	end
	function evolveState!(Ψ, p::Parameters)
		if p.sp.useKrylov
			if p.bhp.isThereDisorder
				propagate!(p, p.pa.tempMatKrylov, Ψ)
			else
				propagate!(p, p.bhp.𝐻, Ψ)
			end
		else
			if p.bhp.isThereDisorder
				Ψ .= normalize(p.pa.tempMatNotKrylov * Ψ)
			else
				Ψ .= normalize(p.pa.expH * Ψ)
			end
		end
	end
	function solveEveryTimeStep(p::Parameters, projectAfterTimeStep)::Array{Array{Complex{Float64},1},1}
		state = copy(p.Ψ₀)
		out = [copy(state)]
		updateTempMatrices!(p)#Generate proper matrices in the memory for disorder etc...
		for i in 2:p.t.steps
			evolveState!(state, p)
			if projectAfterTimeStep
				measurementEffect!(state, p)
			else
				if i in p.t.measIndexes
					measurementEffect!(state, p)
				end
			end
			push!(out, copy(state))
		end
		return out
	end
	function solveLastTimeStep(p::Parameters, projectAfterTimeStep)::Array{Array{Complex{Float64},1},1}
		state = copy(p.Ψ₀)
		updateTempMatrices!(p)
		for i in 2:p.t.steps
			evolveState!(state, p)
			if projectAfterTimeStep
				measurementEffect!(state, p)
			else
				if i in p.t.measIndexes
					measurementEffect!(state, p)
				end
			end
		end
		return [state]
	end
	function MIPT(p::Parameters; onlyLastValue=false, projectAfterTimeStep=false)
		if p.traj > 1000
			out = arrayForEveryThread()
			parameters = parametersForEveryThread(p)
			Threads.@threads for _ in 1:p.traj
				if onlyLastValue
					push!(out[Threads.threadid()], solveLastTimeStep(parameters[Threads.threadid()], projectAfterTimeStep))
				else
					push!(out[Threads.threadid()], solveEveryTimeStep(parameters[Threads.threadid()], projectAfterTimeStep))
				end
			end
			return reduce(vcat, out)
		else
			out = []
			for _ in 1:p.traj
				if onlyLastValue
					push!(out, solveLastTimeStep(p, projectAfterTimeStep))
				else
					push!(out, solveEveryTimeStep(p, projectAfterTimeStep))
				end
			end
			return out
		end
	end
	function arrayForEveryThread()
		a = []
		for _ in 1:Threads.nthreads()
			push!(a, [])
		end
		return a
	end
	function parametersForEveryThread(p::Parameters)
		out = []
		for _ in 1:Threads.nthreads()
			push!(out, deepcopy(p))
		end
		return out
	end
	function savePlotData(x, y, title::String, p::Parameters, initialState::String; notes="")
		path = pwd() * "/Plots/" * title
		mkpath(path)
		io = open(path * "/data.txt", "w")
		println(io, "L = " * string(p.sp.L))
		println(io, "N = " * string(p.sp.N))
		println(io, "cap = " * string(p.sp.cap))
		println(io, "sdim = " * string(p.sdim))
		println(io, "U = " * string(p.bhp.U))
		println(io, "J = " * string(p.bhp.J))
		println(io, "p = " * string(p.pp.p))
		println(io, "f = " * string(p.pp.f))
		println(io, "dt = " * string(p.t.dt))
		println(io, "t = " * string(p.t.times[end]))
		println(io, "Ψ₀ = " * initialState)
		if p.bhp.isThereDisorderInW
			println(io, "w mean = " * string(p.bhp.w))
			println(io, "w stantard deviation = " * string(p.bhp.wσ))
		end
		if p.bhp.isThereDisorderInU
			println(io, "U stantard deviation = " * string(p.bhp.Uσ))
		end
		if p.bhp.isThereDisorderInJ
			println(io, "J stantard deviation = " * string(p.bhp.Jσ))
		end
		println(io, "\n" * notes)
		close(io)
		@save path * "/" * "data.bson" x y
		savefig(path * "/" * title * ".png")
	end
	function expM(M)
		va, U = eigen(M)
		d = U' * M * U
		d .= exp(Diagonal(d))
		return U * d * U'
	end
end
#=
function MIPTOnlyLastValue(p::Parameters)
	out = arrayForEveryThread()
	Threads.@threads for _ in 1:p.traj
		push!(out[Threads.threadid()], solveLastTimeStep(p))
	end
	return reduce(vcat, out)
end
function MIPTProjectAfterEveryTimeStep(p::Parameters)
	out = arrayForEveryThread()
	Threads.@threads for _ in 1:p.traj
		push!(out[Threads.threadid()], solveEveryTimeStepAndProject(p))
	end
	return reduce(vcat, out)
end
function MIPTOnlyLastValueAndProject(p::Parameters)
	out = arrayForEveryThread()
	Threads.@threads for _ in 1:p.traj
		push!(out[Threads.threadid()], solveLastTimeStepAndProject(p))
	end
	return reduce(vcat, out)
end
=#

module MIPTM
	using DifferentialEquations, IterTools, LinearAlgebra, SparseArrays, Plots
	using Distributions, ParametersModule
	using Statistics: mean
	using BSON: @save
	include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl", "OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

	export calcMean, calcMeanAndVar, expVal
	export MIPT
	export singleSubspaceProjectors, onesState, zeroOneState, oneZeroState
	export generateProjectionOperators, generateSingleSite
	export halfBosonNumber, properFluc
	export savePlotData

	function generateProjectionOperators(L, N, cap=N)
		out = []
		for l in 1:L
			oneSiteOperators = []
			for n in 0:min(N, cap)
				push!(oneSiteOperators, projector(L, N, l, n, cap))
			end
			push!(out, oneSiteOperators)
		end
		return out
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
	function singleSubspaceProjectors(L, N, cap=N)
		f(L, N, l) = projector(L, N, l, 1, cap)
		return generateSingleSite(L, N, f)
	end
	function halfBosonNumber(Ψ, L, N, cap=N)
		nₕ = number(L, N, 1, cap)
		for i in 2:Int(round(L/2))
			nₕ .+= number(L, N, i, cap)
		end
		return expVal(Ψ, nₕ^2) - expVal(Ψ, nₕ)^2
	end
	function onesState(L::Int64, cap=L)
		state = zeros(dimensions(L, L, cap))
		state[find_index(ones(Int64, L), cap)] = 1.
		return state
	end
	function oneZeroState(L::Int64, N::Int64, cap=N)
		basisState::Array{Int64,1} = []
		for i in 1:L
			push!(basisState, i%2)
		end
		state = zeros(dimensions(L, N, cap))
		state[find_index(basisState, cap)] = 1.
		return state
	end
	function zeroOneState(L::Int64, N::Int64, cap=N)
		basisState::Array{Int64,1} = [0]
		for i in 1:L-1
			push!(basisState, i%2)
		end
		state = zeros(dimensions(L, N, cap))
		state[find_index(basisState, cap)] = 1.
		return state
	end
	function expVal(s::Array{Complex{Float64},1}, op::SparseMatrixCSC{Float64,Int64})#Jos s on ket
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
		return mean ./ numOfVal
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
        return mean, var
    end
	function properFluc(sol, p)
		nₕ = number(p.sp.L, p.sp.N, 1, p.cap)
		for i in 2:Int(round(p.sp.L/2))
			nₕ .+= number(p.sp.L, p.sp.N, i, p.cap)
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
	function solveEveryTimeStep(p::Parameters, projectAfterTimeStep)
		state = copy(p.sp.Ψ₀)
		out = [state]
		for i in 2:p.t.steps
			if p.bhp.isThereDisorderInW || p.bhp.isThereDisorderInU
				makeDisorderHamiltonian!(p)
				state = propagate(p.bhp.𝐻 .+ p.disorder𝐻, state, p.sdim, p.t.dt)
			else
				state = propagate(p.bhp.𝐻, state, p.sdim, p.t.dt)
			end
			if projectAfterTimeStep
				measurementEffect!(state, p)
			else
				if i in p.t.measIndexes
					measurementEffect!(state, p)
				end
			end
			push!(out, state)
		end
		return out
	end
	function solveLastTimeStep(p::Parameters, projectAfterTimeStep)
		state = copy(p.sp.Ψ₀)
		for i in 2:p.t.steps
			if p.bhp.isThereDisorderInW || p.bhp.isThereDisorderInU
				makeDisorderHamiltonian!(p)
				state .= propagate(p.bhp.𝐻 .+ p.disorder𝐻, state, p.sdim, p.t.dt)
			else
				state .= propagate(p.bhp.𝐻, state, p.sdim, p.t.dt)
			end
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
		out = arrayForEveryThread()
		f = solveEveryTimeStep
		if onlyLastValue
			f = solveLastTimeStep
		end
		Threads.@threads for _ in 1:p.traj
			push!(out[Threads.threadid()], f(p, projectAfterTimeStep))
		end
		return reduce(vcat, out)
	end
	function arrayForEveryThread()
		a = []
		for _ in 1:Threads.nthreads()
			push!(a, [])
		end
		return a
	end
	function savePlotData(x, y, title::String, p::Parameters, initialState::String; notes="")
		path = pwd() * "/Plots/" * title
		mkpath(path)
		io = open(path * "/data.txt", "w")
		println(io, "L = " * string(p.sp.L))
		println(io, "N = " * string(p.sp.N))
		println(io, "cap = " * string(p.cap))
		println(io, "sdim = " * string(p.sdim))
		println(io, "U/J = " * string(p.bhp.U/p.bhp.J))
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
		println(io, "\n" * notes)
		close(io)
		@save path * "/" * "data.bson" x y
		savefig(path * "/" * title * ".png")
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

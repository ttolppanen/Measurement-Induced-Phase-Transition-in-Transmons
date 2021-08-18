module MIPTM
	using DifferentialEquations, IterTools, LinearAlgebra, SparseArrays, Plots
	using Distributions
	using Statistics: mean
	using BSON: @save
	include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl", "OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

	export Parameters, ParametersConstructor, ParametersConstructorWithP, ParametersConstructorWithAny
	export calcMean, calcMeanAndVar, expVal
	export MIPT
	export singleSubspaceProjectors, onesState, zeroOneState, oneZeroState
	export generateProjectionOperators, generateSingleSite
	export halfBosonNumber
	export savePlotData
	export q_next!, q_find_index

	StateType = Union{Array{Float64,1}, Array{Complex{Float64},1}}

	struct TimeData
    	dt::Float64
    	steps::Int64
    	times::Array{Float64, 1}
		measIndexes::Array{Int64, 1}
    	function TimeData(dt::Float64, endTime::Float64, f::Float64)
			times = collect(0.0:dt:endTime)
			steps = length(times)
			measTimes = collect(1/f:1/f:endTime)
			measIndexes = collectMeasIndexes(times, measTimes)
			new(dt, steps, times, measIndexes)
    	end
    end
	function collectMeasIndexes(times, measTimes)
		reversedTimes = reverse(measTimes)
		lastVal = pop!(reversedTimes)
		out = []
		for i in 1:length(times)
			if length(measTimes) == 0
				break
			elseif times[i] > lastVal
				lastVal = pop!(reversedTimes)
				push!(out, i)
			end
		end
		return out
	end
	struct Parameters
		L::Int64 #Number of sites
		N::Int64 #Number of bosons
		cap::Int64 #Max number of bosons on one site
		sdim::Int64
		U::Float64
		J::Float64
		p::Float64 #Probability of measuring one site
		f::Float64 #Frequency of measurements
		t::TimeData #The duration of simulation, and also the time-step
		traj::Int64 #Number of trajectories
		measOp::Array{Any,1}
		𝐻::SparseMatrixCSC{Float64,Int64}
		Ψ₀::Array{Complex{Float64},1}
		μ::Float64 #mean for disorder
		σ::Float64 #stantard deviation
	end
	function ParametersConstructor(;L::Int64, N::Int64, cap=N, dt::Float64, time::Float64,
		traj=1, p::Float64, f::Float64, U::Float64, J::Float64, measOp::Array{Any,1},
		Ψ₀::StateType, sdim::Int64, mean=0.0, stantardDeviation=0.0)
		t = TimeData(dt, time, f)
		HU = interaction(L, N, cap)
		HJ = hopping(L, N, cap)
		d = dimensions(L, N, cap)
		display(d)
		if sdim > d
			display("sdim larger than dimensions! Changed sdim = dimensions.")
			sdim = d
		end
		𝐻 = U .* HU .+ J .* HJ
		return Parameters(L, N, cap, sdim, U, J, p, f, t, traj, measOp, 𝐻, convert(Array{Complex{Float64},1}, Ψ₀), mean, stantardDeviation)
	end
	function ParametersConstructorWithP(p::Parameters, prob::Float64)
		return Parameters(p.L, p.N, p.cap, p.sdim, p.U, p.J, prob, p.f, p.t, p.traj, p.measOp, p.𝐻, p.Ψ₀, p.μ, p.σ)
	end
	function ParametersConstructorWithAny(p::Parameters; L=p.L, cap=p.cap, traj=p.traj)
		return Parameters(L, p.N, cap, p.sdim, p.U, p.J, p.p, p.f, p.t, traj, p.measOp, p.𝐻, p.Ψ₀, p.μ, p.σ)
	end
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
		display(basisState)
		state = zeros(dimensions(L, N, cap))
		state[find_index(basisState, cap)] = 1.
		return state
	end
	function zeroOneState(L::Int64, N::Int64, cap=N)
		basisState::Array{Int64,1} = [0]
		for i in 1:L-1
			push!(basisState, i%2)
		end
		display(basisState)
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
	function measurementEffect!(Ψ, p::Parameters)
		for l in 1:p.L
			if rand(Float64) < p.p  #Does the measurement happen?
				probForProjection = rand(Float64)
				pⱼ = 0 #Probability for a single projection
				for j in 1:length(p.measOp[l])
					if j == length(p.measOp[l])
						projection!(Ψ, p.measOp[l][j])
						break
					else
						pⱼ += projectionProbability(Ψ, p.measOp[l][j])
						if probForProjection < pⱼ
							projection!(Ψ, p.measOp[l][j])
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
	function entanglementAndMeasProbability(p::Parameters, measRates)
		res = []
		halfOfSystems = Int(floor((p.numOfSys / 2)))
		for rate in measRates
			param = NewProbParameters(p=p, Γ=rate)
			@time sol = lastValues(MIPT(param))
			push!(res, calcMean(sol, x -> vonNeumann(x, param.s^halfOfSystems, param.s^(param.numOfSys - halfOfSystems)))[1])
		end
		return res
	end
	function solveEveryTimeStep(p::Parameters, projectAfterTimeStep)
		state = copy(p.Ψ₀)
		out = [state]
		dis = rand(Normal(p.μ, p.σ), p.L)
		#HD = disorder(p.L, p.N, dis=dis)
		#𝐻 = HD .+ p.𝐻
		𝐻 = p.𝐻
		for i in 2:p.t.steps
			state = propagate(𝐻, state, p.sdim, p.t.dt)
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
		state = copy(p.Ψ₀)
		dis = rand(Normal(p.μ, p.σ), p.L)
		#HD = disorder(p.L, p.N, dis=dis)
		𝐻 = p.𝐻
		for i in 2:p.t.steps
			state .= propagate(𝐻, state, p.sdim, p.t.dt)
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
		println(io, "L = " * string(p.L))
		println(io, "N = " * string(p.N))
		println(io, "cap = " * string(p.cap))
		println(io, "sdim = " * string(p.sdim))
		println(io, "U/J = " * string(p.U/p.J))
		println(io, "p = " * string(p.p))
		println(io, "f = " * string(p.f))
		println(io, "dt = " * string(p.t.dt))
		println(io, "t = " * string(p.t.times[end]))
		println(io, "Ψ₀ = " * initialState)
		println(io, "Mean = " * string(p.μ))
		println(io, "Stantard Deviation = " * string(p.σ))
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

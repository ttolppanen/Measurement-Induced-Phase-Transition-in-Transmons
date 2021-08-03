module MIPTM
	using DifferentialEquations, IterTools, LinearAlgebra, SparseArrays, Plots
	using Statistics: mean
	include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl", "OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

	export Parameters, ParametersConstructor, ParametersConstructorWithP
	export calcMean, calcMeanAndVar, expVal
	export MIPT
	export singleSubspaceProjectors, onesState, zeroOneState
	export generateProjectionOperators, generateSingleSite
	export ELVSavefig

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
		sdim::Int64
		p::Float64 #Probability of measuring one site
		f::Float64 #Frequency of measurements
		t::TimeData #The duration of simulation, and also the time-step
		traj::Int64 #Number of trajectories
		measOp::Array{Any,1}
		𝐻::SparseMatrixCSC{Float64,Int64}
		Ψ₀::Array{Complex{Float64},1}
	end
	function ParametersConstructor(;L::Int64, N::Int64, dt::Float64, time::Float64,
		traj=1, p::Float64, f::Float64, U::Float64, J::Float64, measOp::Array{Any,1}, Ψ₀::StateType, sdim::Int64)
		t = TimeData(dt, time, f)
		HU, HJ = split_hamiltonian(L, N)
		𝐻 = U .* HU .+ J .* HJ
		return Parameters(L, N, sdim, p, f, t, traj, measOp, 𝐻, convert(Array{Complex{Float64},1}, Ψ₀))
	end
	function ParametersConstructorWithP(p::Parameters, prob::Float64)
		return Parameters(p.L, p.N, p.sdim, prob, p.f, p.t, p.traj, p.measOp, p.𝐻, p.Ψ₀)
	end
	function generateProjectionOperators(L, N)
		out = []
		for l in 1:L
			oneSiteOperators = []
			for n in 0:N
				push!(oneSiteOperators, projector(L, N, l, n))
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
	function singleSubspaceProjectors(L, N)
		f(L, N, l) = projector(L, N, l, 1)
		return generateSingleSite(L, N, f)
	end
	function onesState(L::Int64)
		state = zeros(dimension(L, L))
		state[find_index(ones(Int64, L))] = 1.;
		return state
	end
	function zeroOneState(L::Int64, N::Int64)
		basisState::Array{Int64,1} = []
		for i in 1:L
			push!(basisState, i%2)
		end
		state = zeros(dimension(L, N))
		state[find_index(basisState)] = 1.;
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
	function solveEveryTimeStep(p::Parameters)
		state = copy(p.Ψ₀)
		out = [state]
		for i in 2:p.t.steps
			state = propagate(p.𝐻, state, p.sdim, p.t.dt)
			if i in p.t.measIndexes
				measurementEffect!(state, p)
			end
			push!(out, state)
		end
		return out
	end
	function solveEveryTimeStepAndProject(p::Parameters) #Projection after every time-step
		state = copy(p.Ψ₀)
		out = [state]
		for i in 2:p.t.steps
			state = propagate(p.𝐻, state, p.sdim, p.t.dt)
			measurementEffect!(state, p)
			push!(out, state)
		end
		return out
	end
	function solveLastTimeStep(p::Parameters)
		state = copy(p.Ψ₀)
		for i in 2:p.t.steps
			state .= propagate(p.𝐻, state, p.sdim, p.t.dt)
			if i in p.t.measIndexes
				measurementEffect!(state, p)
			end
		end
		return [state]
	end
	function solveLastTimeStepAndProject(p::Parameters)
		state = copy(p.Ψ₀)
		for i in 2:p.t.steps
			state .= propagate(p.𝐻, state, p.sdim, p.t.dt)
			measurementEffect!(state, p)
		end
		return [state]
	end
	function MIPT(p::Parameters; onlyLastValue=false, projectAfterTimeStep=false)
		out = arrayForEveryThread()

		f = solveEveryTimeStep
		if onlyLastValue && projectAfterTimeStep
			f = solveLastTimeStepAndProject
		elseif onlyLastValue
			f = solveLastTimeStep
		elseif projectAfterTimeStep
			f = solveEveryTimeStepAndProject
		end

		Threads.@threads for _ in 1:p.traj
			push!(out[Threads.threadid()], f(p))
		end
		return reduce(vcat, out)
	end
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
	function arrayForEveryThread()
		a = []
		for _ in 1:Threads.nthreads()
			push!(a, [])
		end
		return a
	end
	function ELVSavefig(;title::String, p::Parameters, initialState::String, notes="")
		folderName = "ELV_" * title
		path = pwd() * "/Plots/" * folderName
		mkpath(path)
		io = open(path * "/data.txt", "w")
		println(io, "L = " * string(p.L))
		println(io, "N = " * string(p.N))
		println(io, "sdim = " * string(p.sdim))
		println(io, "p = " * string(p.p))
		println(io, "f = " * string(p.f))
		println(io, "dt = " * string(p.t.dt))
		println(io, "t = " * string(p.t.times[end]))
		println(io, "Ψ₀ = " * initialState)
		println(io, "\n" * notes)
		close(io)
		savefig(pwd() * "/Plots/" * folderName * "/" * folderName * ".png")
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

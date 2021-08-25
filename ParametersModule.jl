module ParametersModule

	using SparseArrays, Distributions
	include.(["OllisCode/Operators.jl", "OllisCode/Basis.jl"])

	export SystemParameters, ProjectionParameters, BoseHubbardParameters, Parameters
	export makeDisorderHamiltonian!, shouldUseKrylov

	StateType = Union{Array{Float64,1}, Array{Complex{Float64},1}}
	MatrixType = Union{Array{Float64,2}, Array{Complex{Float64},2}, SparseMatrixCSC{Float64,Int64}}

	struct TimeData
		dt::Float64
		steps::Int64
		times::Array{Float64, 1}
		measIndexes::Array{Int64, 1}
		function TimeData(dt::Float64, endTime::Float64, f::Float64)
			times = collect(0.0:dt:endTime)
			steps = length(times)
			measTimes = collect(1/f:1/f:endTime)
			measIndexes = collectMeasIndexes(times, measTimes) #I think these are the indexes when measurements should occur, so that dt is always constant...
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
	mutable struct SystemParameters
		L::Int64 #Number of sites
		N::Int64 #Number of bosons
		cap::Int64
		dim::Int64
		useKrylov::Bool
		function SystemParameters(;L::Int64, N::Int64, cap=N)
			dim = dimensions(L, N, cap=cap)
			useKrylov = shouldUseKrylov(dim)
			new(L, N, cap, dim, useKrylov)
		end
	end
	mutable struct ProjectionParameters
		p::Float64 #Probability of measuring one site
		f::Float64 #Frequency of measurements
		projOp::Array{Any,1} #Projection operators
		function ProjectionParameters(;p::Float64, f::Float64, projOp::Array{Any,1})
			new(p, f, projOp)
		end
	end
	mutable struct BoseHubbardParameters
		w::Float64
		wσ::Float64 #For disorder
		isThereDisorderInW::Bool
		U::Float64
		Uσ::Float64 #For disorder
		isThereDisorderInU::Bool
		J::Float64
		𝐻::MatrixType
		function BoseHubbardParameters(;L::Int64, N::Int64, cap=N,
				w::Float64=0.0, wσ::Float64=0.0,
				U::Float64, Uσ::Float64=0.0, J::Float64=1.0)

			isThereDisorderInW = wσ != 0.0
			isThereDisorderInU = Uσ != 0.0
			HJ = hopping(L, N, cap=cap)
			𝐻 = J .* HJ
			if !isThereDisorderInU
				𝐻 .+= U * interaction(L, N, cap=cap)
			end
			new(w, wσ, isThereDisorderInW, U, Uσ, isThereDisorderInU, J, 𝐻)
		end
	end
	mutable struct Parameters
		sp::SystemParameters
		pp::ProjectionParameters
		bhp::BoseHubbardParameters
		sdim::Int64
		t::TimeData #The duration of simulation, and also the time-step
		traj::Int64 #Number of trajectories
		disorder𝐻::MatrixType #temp matrix for disorder
		Ψ₀::Array{Complex{Float64},1}
		function Parameters(;sp::SystemParameters, pp::ProjectionParameters,
					bhp::BoseHubbardParameters, dt::Float64, time::Float64,
					traj=1, sdim::Int64, Ψ₀::StateType)

			t = TimeData(dt, time, pp.f)
			display(sp.dim)
			if sdim > sp.dim
				display("sdim larger than dimensions! Changed sdim = dimensions.")
				sdim = sp.dim
			end
			if sp.useKrylov
				disorder𝐻 = spzeros(sp.dim, sp.dim)
			else
				disorder𝐻 = zeros(sp.dim, sp.dim)
			end
			new(sp, pp, bhp, sdim, t, traj, disorder𝐻, convert(Array{Complex{Float64},1}, Ψ₀))
		end
	end
	function makeDisorderHamiltonian!(p::Parameters)
		if p.bhp.isThereDisorderInW && p.bhp.isThereDisorderInU
			p.disorder𝐻 .= disorder(p.sp.L, p.sp.N, p.sp.cap, dis = disorderForW(p.bhp, p.sp.L))  .+ interaction(p.sp.L, p.sp.N, p.sp.cap, dis = disorderForU(p.bhp, p.sp.L))
		elseif p.bhp.isThereDisorderInW
			p.disorder𝐻 .= disorder(p.sp.L, p.sp.N, p.sp.cap, dis = disorderForW(p.bhp, p.sp.L))
		elseif p.bhp.isThereDisorderInU
			p.disorder𝐻 .= interaction(p.sp.L, p.sp.N, p.sp.cap, dis = disorderForU(p.bhp, p.sp.L))
		else
			error("Trying to create a Hamiltonian with disorder when there isn't supposed to be any...")
		end
	end
	#=
	function ParametersConstructorWithP(p::Parameters, prob::Float64)
		return Parameters(p.L, p.N, p.cap, p.sdim, p.U, p.J, prob, p.f, p.t, p.traj, p.measOp, p.𝐻, p.Ψ₀, p.μ, p.σ)
	end
	function ParametersConstructorWithAny(p::Parameters; L=p.L, cap=p.cap, traj=p.traj)
		return Parameters(L, p.N, cap, p.sdim, p.U, p.J, p.p, p.f, p.t, traj, p.measOp, p.𝐻, p.Ψ₀, p.μ, p.σ)
	end
	=#

	function returnDisorder(L::Int64, μ, σ)
		return rand(Normal(μ, σ), L)
	end
	function disorderForW(bhp::BoseHubbardParameters, L)
		returnDisorder(L, bhp.w, bhp.wσ)
	end
	function disorderForU(bhp::BoseHubbardParameters, L)
		returnDisorder(L, bhp.U, bhp.Uσ)
	end
	function shouldUseKrylov(dim)
		#return dim >= 2002
		return false
	end
end

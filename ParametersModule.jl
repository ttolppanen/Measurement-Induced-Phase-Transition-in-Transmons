module ParametersModule

	using SparseArrays, Distributions, LinearAlgebra, InteractiveUtils
	include.(["OllisCode/Operators.jl", "OllisCode/Basis.jl"])

	export SystemParameters, ProjectionParameters, BoseHubbardParameters, Parameters
	export updateTempMatrices!, shouldUseKrylov

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
			if cap > N
				display("cap > N, setting cap = N")
				cap = N
			end
			dim = dimensions(L, N, cap=cap)
			useKrylov = dim > 1000
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
	struct BoseHubbardParameters
		w::Float64
		wσ::Float64 #For disorder
		isThereDisorderInW::Bool
		U::Float64
		Uσ::Float64 #For disorder
		isThereDisorderInU::Bool
		J::Float64
		Jσ::Float64
		isThereDisorderInJ::Bool
		isThereDisorder::Bool
		𝐻::SparseMatrixCSC{Float64,Int64}
		function BoseHubbardParameters(;sp::SystemParameters,
				w::Float64=0.0, wσ::Float64=0.0,
				U::Float64=0.0, Uσ::Float64=0.0,
				J::Float64=1.0, Jσ::Float64=0.0)

			isThereDisorderInW = wσ != 0.0
			isThereDisorderInU = Uσ != 0.0
			isThereDisorderInJ = Jσ != 0.0
			isThereDisorder = isThereDisorderInW || isThereDisorderInU || isThereDisorderInJ
			𝐻 = spzeros(sp.dim, sp.dim)
			if !isThereDisorderInJ
				𝐻 .+= J .* hopping(sp.L, sp.N, sp.dim, cap=sp.cap)
			end
			if !isThereDisorderInU
				𝐻 .+= U .* interaction(sp.L, sp.N, sp.dim, cap=sp.cap)
			end
			new(w, wσ, isThereDisorderInW, U, Uσ, isThereDisorderInU, J, Jσ, isThereDisorderInJ, isThereDisorder, 𝐻)
		end
	end
	mutable struct PreAllocated
		tempMatKrylov::Union{SparseMatrixCSC{Float64,Int64}, Nothing} #Matrices for disordered hamiltonian or the time evolution...
		tempMatNotKrylov::Union{Array{Complex{Float64},2}, Nothing}
		expH::Union{Array{Complex{Float64},2}, Nothing}
		V::Union{Array{Complex{Float64},2}, Nothing} #for krylov!
		h::Union{Array{Complex{Float64},2}, Nothing}
		w::Union{Array{Complex{Float64},1}, Nothing}
		function PreAllocated(dim, sdim, dt, 𝐻, useKrylov)
			if useKrylov
				tempMatKrylov = spzeros(dim, dim)
				V = zeros(ComplexF64, dim, sdim)
				h = zeros(ComplexF64, sdim, sdim)
				w = zeros(ComplexF64, dim)
				new(tempMatKrylov, nothing, nothing, V, h, w)
			else
				tempMatNotKrylov = zeros(dim, dim)
				expH = expM(-1im .* dt .* Matrix(𝐻))
				new(nothing, tempMatNotKrylov, expH, nothing, nothing, nothing)
			end
		end
	end
	mutable struct Parameters
		sp::SystemParameters
		pp::ProjectionParameters
		bhp::BoseHubbardParameters
		pa::PreAllocated
		sdim::Int64
		t::TimeData #The duration of simulation, and also the time-step
		traj::Int64 #Number of trajectories
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
			pa = PreAllocated(sp.dim, sdim, dt, bhp.𝐻, sp.useKrylov)
			new(sp, pp, bhp, pa, sdim, t, traj, convert(Array{Complex{Float64},1}, Ψ₀))
		end
	end
	function updateTempMatrices!(p::Parameters)
		if p.bhp.isThereDisorder
			if p.sp.useKrylov
				makeDisorderHamiltonian!(p::Parameters)
			else
				makeDisorderedTimeEvolution!(p::Parameters)
			end
		end
	end
	function makeDisorderHamiltonian!(p::Parameters)
		p.pa.tempMatKrylov .= copy(p.bhp.𝐻)
		if p.bhp.isThereDisorderInW
			p.pa.tempMatKrylov .+= disorder(p.sp.L, p.sp.N, p.sp.dim, cap=p.sp.cap, dis = disorderForW(p.bhp, p.sp.L))
		end
		if p.bhp.isThereDisorderInU
			p.pa.tempMatKrylov .+= interaction(p.sp.L, p.sp.N, p.sp.dim, cap=p.sp.cap, dis = disorderForU(p.bhp, p.sp.L))
		end
		if p.bhp.isThereDisorderInJ
			p.pa.tempMatKrylov .+= hopping(p.sp.L, p.sp.N, p.sp.dim, cap=p.sp.cap, dis = disorderForJ(p.bhp, p.sp.L))
		end
	end
	function makeDisorderedTimeEvolution!(p::Parameters)
		mat = copy(p.bhp.𝐻)
		if p.bhp.isThereDisorderInW
			mat .+= disorder(p.sp.L, p.sp.N, p.sp.dim, cap=p.sp.cap, dis = disorderForW(p.bhp, p.sp.L))
		end
		if p.bhp.isThereDisorderInU
			mat .+= interaction(p.sp.L, p.sp.N, p.sp.dim, cap=p.sp.cap, dis = disorderForU(p.bhp, p.sp.L))
		end
		if p.bhp.isThereDisorderInJ
			mat .+= hopping(p.sp.L, p.sp.N, p.sp.dim, cap=p.sp.cap, dis = disorderForJ(p.bhp, p.sp.L))
		end
		p.pa.tempMatNotKrylov .= expM(-1im .* p.t.dt .* Matrix(mat))
		nothing
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
		return returnDisorder(L, bhp.w, bhp.wσ)
	end
	function disorderForU(bhp::BoseHubbardParameters, L)
		return returnDisorder(L, bhp.U, bhp.Uσ)
	end
	function disorderForJ(bhp::BoseHubbardParameters, L)
		return returnDisorder(L, bhp.J, bhp.Jσ)
	end
	function expM(M::Array{Complex{Float64},2})::Array{Complex{Float64},2}
		va, U = eigen(M)
		d = U' * M * U
		d .= exp(Diagonal(d))
		return U * d * U'
	end
end

module MIPTM
	using DifferentialEquations, IterTools, LinearAlgebra
	using Statistics: mean

	export Parameters, kronForMany, calcMean, ensSolToList, expVal

	struct TimeData
    	dt::Float64
    	Δt::Tuple{Float64,Float64}
    	times::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    	function TimeData(startTime::Float64, dt::Float64, endTime::Float64)
    		new(dt, (startTime, endTime), startTime:dt:endTime)
    	end
    end
	struct Operators
		n::Array{Complex{Float64},2}
		nAll::Array{Complex{Float64},2}
		n𝐼::Array{Array{Complex{Float64},2},1}
		a::Array{Complex{Float64},2}
		ad::Array{Complex{Float64},2}
		𝐼::Array{Complex{Float64},2}
		function Operators(s::Int64, numOfSys::Int64)
			𝐼 = makeI(s)
		    a = lowOp(s)
		    ad = copy(a')
		    n = ad*a
			n𝐼 = listOfOperators(n, numOfSys, 𝐼)
		    nAll =  sum(n𝐼)
			new(n, nAll, n𝐼, a, ad, 𝐼)
		end
	end
	struct Parameters
		numOfSys::Int64 #Number of systems
		s::Int64 #Max size of one system
		dim::Int64 #The dimension of the total system
		p::Float64 #Probability of measuring one site
		f::Float64 #Frequency of these measurements
		t::TimeData #The duration of simulation, and also the time-step
		traj::Int64 #Number of trajectories
		atol::Float64
		rtol::Float64
		op::Operators #Some usefull operators
		𝐻::Array{Complex{Float64},2}
		function Parameters(;numOfSys::Int64,s::Int64,t::Tuple{Float64,Float64,Float64},
			traj=1, atol=1e-3, rtol=1e-3, p::Float64, f::Float64, ω::Float64, U::Float64, J::Float64)
			dim = s^numOfSys
			op = Operators(s, numOfSys)
			t = TimeData(t[1], t[2], t[3])
			new(numOfSys, s, dim, p, f, t, traj, atol, rtol, op, boseHubbard(ω=ω, U=U, J=J, n=op.n, a=op.a, 𝐼=op.𝐼, numOfSys=numOfSys))
		end
	end
	function lowOp(s::Int64) #â
		a = zeros(s, s)
		for i in 1:s-1
			a[i, i + 1] = sqrt(i)
		end
		complex(a)
	end
	function makeI(size::Int64)
        𝐼 = (1.0 + 0.0*im)*Matrix(I, size, size)
    end
	function boseHubbard(;ω::Float64, U::Float64, J::Float64, n::Array{Complex{Float64},2}, a::Array{Complex{Float64},2}, 𝐼::Array{Complex{Float64},2}, numOfSys::Int64)
		nᵢ = kronForMany(n, 𝐼, 1, numOfSys)
		𝐼All = kronForMany(𝐼, 𝐼, 1, numOfSys)
		H = ω*nᵢ - U*0.5*nᵢ*(nᵢ-𝐼All)
		for i in 2:numOfSys
			nᵢ = kronForMany(n, 𝐼, i, numOfSys)
			aᵢ₋₁ = kronForMany(a, 𝐼, i - 1, numOfSys)
			aᵢ = kronForMany(a, 𝐼, i, numOfSys)
			H .+= ω*nᵢ - U*0.5*nᵢ*(nᵢ-𝐼All) + J*(aᵢ₋₁*aᵢ' + aᵢ₋₁'*aᵢ)
		end
		H
	end
	function partialTrace(ρ::Array{Complex{Float64},2}, aDim::Int64, bDim::Int64; traceOverB::Bool=true)::Array{Complex{Float64},2}
        if traceOverB
            A = complex(zeros(aDim, aDim))
            for i in 1:aDim
                iᵨ = 1 + (i - 1)*bDim
                for j in 1:aDim
                    jᵨ = 1 + (j - 1)*bDim
                    for d in 0:bDim-1
                        A[i, j] += ρ[iᵨ+d, jᵨ+d]
                    end
                end
            end
            A
        else
            B = complex(zeros(bDim, bDim))
            for i in 1:bDim
                for j in 1:bDim
                    for d in 0:aDim-1
                        B[i, j] += ρ[i + d*bDim, j + d*bDim]
                    end
                end
            end
            B
        end
    end
	function kronForMany(m::Array{Complex{Float64},2}, 𝐼, index, numOfSys)::Array{Complex{Float64},2}
        if index == numOfSys
            s = m
        else
            s = 𝐼
        end
        for i in reverse(1:numOfSys-1)
            if i == index
                s = kron(m, s)
            else
                s = kron(𝐼, s)
            end
        end
        s
    end
	function kronForMany(m::Union{Array{Array{Complex{Float64},2},1}, Array{Array{Complex{Float64},1},1}})
        s = m[end]
        for (isFirst, mᵢ) in flagfirst(reverse(m))
            if isFirst
            else
                s = kron(mᵢ, s)
            end
        end
        s
    end
	function listOfOperators(op::Array{Complex{Float64},2}, numOfSys::Int64, 𝐼::Array{Complex{Float64},2})::Array{Array{Complex{Float64},2},1}
        res = []
        for i in 1:numOfSys
            push!(res, kronForMany(op, 𝐼, i, numOfSys))
        end
        res
    end
	function expVal(s::Array{Complex{Float64},1}, op::Array{Complex{Float64},2})#Jos s on ket
		real(s' * op * s)
	end
	function expVal(ρ::Array{Complex{Float64},2}, op::Array{Complex{Float64},2})#Tiheysoperaattorille
		real(tr(op*ρ))
	end
	function expVal(ρ::Array{Complex{Float64},2}, op::Array{Complex{Float64},2}, mPA1::Array{Complex{Float64},2})#Tiheysoperaattorille
		mul!(mPA1, op, ρ)
		real(tr(mPA1))
	end
	function calcMean(sol, f::Function)
		numOfVal = length(sol)
		mean = f.(sol[1])
        for i in 2:numOfVal
            mean .+= f.(sol[i])
        end
		mean./numOfVal
    end
	function calcMeanAndVar(sol, f::Function)
		numOfVal = length(sol)
		fVal = f.(sol[1])
		mean = fVal
		var = fVal.^2
        for i in 2:numOfVal
            fVal = f.(sol[i])
            mean .+= fVal
            var .+= fVal.^2
        end
		mean .= mean./numOfVal
		var .= var./numOfVal .- mean.^2
        mean, var
    end
	function ensSolToList(ensSol::EnsembleSolution)::Array{Array{Array{Complex{Float64},1},1},1}
        res = []
		for sol in ensSol
            push!(res, [i for i in sol.u])
        end
		res
    end
end

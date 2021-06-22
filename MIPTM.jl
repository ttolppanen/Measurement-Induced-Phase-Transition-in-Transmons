module MIPTM
	using DifferentialEquations, IterTools, LinearAlgebra
	using Statistics: mean

	export Parameters, ParametersConstructor, NewProbParameters, kronForMany, calcMean, calcMeanAndVar, ensSolToList, expVal, schrodinger
	export MIPT, vonNeumann, entanglementAndMeasProbability

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
		measOp::Array{Array{Array{Complex{Float64},2},1},1}
		a::Array{Complex{Float64},2}
		ad::Array{Complex{Float64},2}
		𝐼::Array{Complex{Float64},2}
		function Operators(s::Int64, numOfSys::Int64, measOp::Array{Array{Complex{Float64},2},1})
			𝐼 = makeI(s)
		    a = lowOp(s)
		    ad = copy(a')
		    n = ad*a
			n𝐼 = listOfOperators(n, numOfSys, 𝐼)
			op = makeSetOfMeasurementOperators(measOp, numOfSys, 𝐼)
		    nAll =  sum(n𝐼)
			new(n, nAll, n𝐼, op, a, ad, 𝐼)
		end
	end
	struct Parameters
		numOfSys::Int64 #Number of systems
		s::Int64 #Max size of one system
		dim::Int64 #The dimension of the total system
		Γ::Float64 #Probability of measuring one site
		t::TimeData #The duration of simulation, and also the time-step
		traj::Int64 #Number of trajectories
		atol::Float64
		rtol::Float64
		op::Operators #Some usefull operators
		𝐻::Array{Complex{Float64},2}
		Ψ₀::Array{Complex{Float64},1}
	end
	function ParametersConstructor(;t::Tuple{Float64,Float64,Float64},
		traj::Int64, atol=1e-3, rtol=1e-3, Γ::Float64, ω::Float64,
		U::Float64, J::Float64, Ψ₀::Array{Array{Complex{Float64},1},1},
		measOp::Array{Array{Complex{Float64},2},1})
		numOfSys = length(Ψ₀)
		s = length(Ψ₀[1])
		dim = s^numOfSys
		op = Operators(s, numOfSys, measOp)
		t = TimeData(t[1], t[2], t[3])
		Parameters(numOfSys, s, dim, Γ, t, traj, atol, rtol, op, boseHubbard(ω=ω, U=U, J=J, n=op.n, a=op.a, 𝐼=op.𝐼, numOfSys=numOfSys), kronForMany(Ψ₀))
	end
	function NewProbParameters(;p::Parameters, Γ::Float64)
		Parameters(p.numOfSys, p.s, p.dim, Γ, p.t, p.traj, p.atol, p.rtol, p.op, p.𝐻, p.Ψ₀)
	end
	function makeSetOfMeasurementOperators(operators, numOfSys, 𝐼)
		measOp = []
		for i in 1:numOfSys
			oneSiteOperators = []
			for j in 1:length(operators)
				push!(oneSiteOperators, kronForMany(operators[j], 𝐼, i, numOfSys))
			end
			push!(measOp, oneSiteOperators)
		end
		measOp
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
	function vonNeumann(Ψ::Array{Complex{Float64},1}, aDim::Int64, bDim::Int64)
		ρₐ = partialTrace(Ψ*Ψ', aDim, bDim)
		F = svd(ρₐ)
		-dot(real(F.S), vonNeumannlog.(real(F.S)))
	end
	function vonNeumannlog(x)
		if x == 0
			return 1
		else
			return log(x)
		end
	end
	function ensSolToList(ensSol::EnsembleSolution, times)::Array{Array{Array{Complex{Float64},1},1},1}
        res = []
		for sol in ensSol
			oneSol = []
			for t in times
				push!(oneSol, sol(t))
			end
			push!(res, oneSol)
        end
		res
    end
	function schrodinger(Ψ, p, t)
		-1im * p.𝐻 * Ψ
	end
	function measurementEffect(integrator)
		Ψ = integrator.u
		p = integrator.p
		dt = integrator.t - integrator.tprev
		for i in 1:p.numOfSys
			if rand(Float64) < p.Γ*dt  #Does the measurement happen?
				probForProjection = rand(Float64)
				pⱼ = 0 #Probability for a single projection
				for j in 1:length(p.op.measOp[i])
					if j == length(p.op.measOp[i])
						projection!(Ψ, p.op.measOp[i][j])
						break
					else
						pⱼ += projectionProbability(Ψ, p.op.measOp[i][j])
						if probForProjection < pⱼ
							projection!(Ψ, p.op.measOp[i][j])
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
		expVal(Ψ, op' * op)
	end
	function MIPT(p::Parameters)
		condition(u, t, integrator) = true
		cb = DiscreteCallback(condition, measurementEffect, save_positions=(true,true))
		prob = ODEProblem(schrodinger, p.Ψ₀, p.t.Δt, p, saveat = p.t.dt)
		enProb = EnsembleProblem(prob, safetycopy=true)
		sol = solve(enProb, SRA1(), abstol=p.atol, rtol=p.rtol, EnsembleThreads(), trajectories=p.traj, callback=cb)
		ensSolToList(sol, p.t.times)
	end
	function lastValues(sol)
		res = []
		for i in sol
			push!(res, [last(i)])
		end
		res
	end
	function entanglementAndMeasProbability(p::Parameters, measRates)
		res = []
		halfOfSystems = Int(floor((p.numOfSys / 2)))
		for rate in measRates
			param = NewProbParameters(p=p, Γ=rate)
			@time sol = lastValues(MIPT(param))
			push!(res, calcMean(sol, x -> vonNeumann(x, param.s^halfOfSystems, param.s^(param.numOfSys - halfOfSystems)))[1])
		end
		res
	end
end

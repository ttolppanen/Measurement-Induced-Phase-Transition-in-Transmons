module MIPTMLegacy
	using DifferentialEquations, IterTools, LinearAlgebra
	using Statistics: mean

	export Parameters, ParametersConstructor, NewProbParameters, kronForMany, calcMean, calcMeanAndVar, ensSolToList, expVal, schrodinger
	export MIPT, vonNeumann, entanglementAndMeasProbability, vonNeumannHalfOfSystem, lowOp
	export makeXProjectors

	struct TimeData
    	dt::Float64
    	Īt::Tuple{Float64,Float64}
    	times::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    	function TimeData(startTime::Float64, dt::Float64, endTime::Float64)
			new(dt, (startTime, endTime), startTime:dt:endTime)
    	end
    end
	struct Operators
		n::Array{Complex{Float64},2}
		nAll::Array{Complex{Float64},2}
		nš¼::Array{Array{Complex{Float64},2},1}
		measOp::Array{Array{Array{Complex{Float64},2},1},1}
		a::Array{Complex{Float64},2}
		ad::Array{Complex{Float64},2}
		š¼::Array{Complex{Float64},2}
		function Operators(s::Int64, numOfSys::Int64, measOp::Array{Array{Complex{Float64},2},1})
			š¼ = makeI(s)
		    a = lowOp(s)
		    ad = copy(a')
		    n = ad*a
			nš¼ = listOfOperators(n, numOfSys, š¼)
			op = makeSetOfMeasurementOperators(measOp, numOfSys, š¼)
		    nAll =  sum(nš¼)
			new(n, nAll, nš¼, op, a, ad, š¼)
		end
	end
	struct Parameters
		numOfSys::Int64 #Number of systems
		s::Int64 #Max size of one system
		dim::Int64 #The dimension of the total system
		Ī::Float64 #Probability of measuring one site
		t::TimeData #The duration of simulation, and also the time-step
		traj::Int64 #Number of trajectories
		atol::Float64
		rtol::Float64
		op::Operators #Some usefull operators
		š»::Array{Complex{Float64},2}
		ĪØā::Array{Complex{Float64},1}
	end
	function ParametersConstructor(;t::Tuple{Float64,Float64,Float64},
		traj::Int64, atol=1e-3, rtol=1e-3, Ī::Float64, Ļ::Float64,
		U::Float64, J::Float64, ĪØā::Array{Array{Complex{Float64},1},1},
		measOp::Array{Array{Complex{Float64},2},1})
		numOfSys = length(ĪØā)
		s = length(ĪØā[1])
		dim = s^numOfSys
		op = Operators(s, numOfSys, measOp)
		t = TimeData(t[1], t[2], t[3])
		Parameters(numOfSys, s, dim, Ī, t, traj, atol, rtol, op, boseHubbard(Ļ=Ļ, U=U, J=J, n=op.n, a=op.a, š¼=op.š¼, numOfSys=numOfSys), kronForMany(ĪØā))
	end
	function NewProbParameters(;p::Parameters, Ī::Float64)
		Parameters(p.numOfSys, p.s, p.dim, Ī, p.t, p.traj, p.atol, p.rtol, p.op, p.š», p.ĪØā)
	end
	function makeSetOfMeasurementOperators(operators, numOfSys, š¼)
		measOp = []
		for i in 1:numOfSys
			oneSiteOperators = []
			for j in 1:length(operators)
				push!(oneSiteOperators, kronForMany(operators[j], š¼, i, numOfSys))
			end
			push!(measOp, oneSiteOperators)
		end
		measOp
	end
	function lowOp(s::Int64) #aĢ
		a = zeros(s, s)
		for i in 1:s-1
			a[i, i + 1] = sqrt(i)
		end
		complex(a)
	end
	function makeI(size::Int64)
        š¼ = (1.0 + 0.0*im)*Matrix(I, size, size)
    end
	function boseHubbard(;Ļ::Float64, U::Float64, J::Float64, n::Array{Complex{Float64},2}, a::Array{Complex{Float64},2}, š¼::Array{Complex{Float64},2}, numOfSys::Int64)
		nįµ¢ = kronForMany(n, š¼, 1, numOfSys)
		š¼All = kronForMany(š¼, š¼, 1, numOfSys)
		H = Ļ*nįµ¢ - U*0.5*nįµ¢*(nįµ¢-š¼All)
		for i in 2:numOfSys
			nįµ¢ = kronForMany(n, š¼, i, numOfSys)
			aįµ¢āā = kronForMany(a, š¼, i - 1, numOfSys)
			aįµ¢ = kronForMany(a, š¼, i, numOfSys)
			H .+= Ļ*nįµ¢ - U*0.5*nįµ¢*(nįµ¢-š¼All) + J*(aįµ¢āā*aįµ¢' + aįµ¢āā'*aįµ¢)
		end
		H
	end
	function partialTrace(Ļ::Array{Complex{Float64},2}, aDim::Int64, bDim::Int64; traceOverB::Bool=true)::Array{Complex{Float64},2}
        if traceOverB
            A = complex(zeros(aDim, aDim))
            for i in 1:aDim
                iįµØ = 1 + (i - 1)*bDim
                for j in 1:aDim
                    jįµØ = 1 + (j - 1)*bDim
                    for d in 0:bDim-1
                        A[i, j] += Ļ[iįµØ+d, jįµØ+d]
                    end
                end
            end
            A
        else
            B = complex(zeros(bDim, bDim))
            for i in 1:bDim
                for j in 1:bDim
                    for d in 0:aDim-1
                        B[i, j] += Ļ[i + d*bDim, j + d*bDim]
                    end
                end
            end
            B
        end
    end
	function kronForMany(m::Array{Complex{Float64},2}, š¼, index, numOfSys)::Array{Complex{Float64},2}
        if index == numOfSys
            s = m
        else
            s = š¼
        end
        for i in reverse(1:numOfSys-1)
            if i == index
                s = kron(m, s)
            else
                s = kron(š¼, s)
            end
        end
        s
    end
	function kronForMany(m::Union{Array{Array{Complex{Float64},2},1}, Array{Array{Complex{Float64},1},1}})
        s = m[end]
        for (isFirst, mįµ¢) in flagfirst(reverse(m))
            if isFirst
            else
                s = kron(mįµ¢, s)
            end
        end
        s
    end
	function listOfOperators(op::Array{Complex{Float64},2}, numOfSys::Int64, š¼::Array{Complex{Float64},2})::Array{Array{Complex{Float64},2},1}
        res = []
        for i in 1:numOfSys
            push!(res, kronForMany(op, š¼, i, numOfSys))
        end
        res
    end
	function expVal(s::Array{Complex{Float64},1}, op::Array{Complex{Float64},2})#Jos s on ket
		real(s' * op * s)
	end
	function expVal(Ļ::Array{Complex{Float64},2}, op::Array{Complex{Float64},2})#Tiheysoperaattorille
		real(tr(op*Ļ))
	end
	function expVal(Ļ::Array{Complex{Float64},2}, op::Array{Complex{Float64},2}, mPA1::Array{Complex{Float64},2})#Tiheysoperaattorille
		mul!(mPA1, op, Ļ)
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
	function vonNeumann(ĪØ::Array{Complex{Float64},1}, aDim::Int64, bDim::Int64)
		Ļā = partialTrace(ĪØ*ĪØ', aDim, bDim)
		F = svd(Ļā)
		-dot(real(F.S), vonNeumannlog.(real(F.S)))
	end
	function vonNeumannlog(x)
		if x == 0
			return 1
		else
			return log(x)
		end
	end
	function vonNeumannHalfOfSystem(ĪØ::Array{Complex{Float64},1}, p::Parameters)
		halfOfSystems = Int(floor((p.numOfSys / 2)))
		vonNeumann(ĪØ, p.s^halfOfSystems, p.s^(p.numOfSys - halfOfSystems))
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
	function schrodinger(ĪØ, p, t)
		-1im * p.š» * ĪØ
	end
	function measurementEffect(integrator)
		ĪØ = integrator.u
		p = integrator.p
		dt = integrator.t - integrator.tprev
		for i in 1:p.numOfSys
			if rand(Float64) < p.Ī  #Does the measurement happen?
				probForProjection = rand(Float64)
				pā±¼ = 0 #Probability for a single projection
				for j in 1:length(p.op.measOp[i])
					if j == length(p.op.measOp[i])
						projection!(ĪØ, p.op.measOp[i][j])
						break
					else
						pā±¼ += projectionProbability(ĪØ, p.op.measOp[i][j])
						if probForProjection < pā±¼
							projection!(ĪØ, p.op.measOp[i][j])
							break
						end
					end
				end
			end
		end
	end
	function projection!(ĪØ, op)
		ĪØ .= op * ĪØ
		ĪØ ./= norm(ĪØ)
	end
	function projectionProbability(ĪØ, op)
		expVal(ĪØ, op' * op)
	end
	function MIPT(p::Parameters)
		condition(u, t, integrator) = true
		cb = DiscreteCallback(condition, measurementEffect, save_positions=(true,true))
		prob = ODEProblem(schrodinger, p.ĪØā, p.t.Īt, p, saveat = p.t.dt)
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
			param = NewProbParameters(p=p, Ī=rate)
			@time sol = lastValues(MIPT(param))
			push!(res, calcMean(sol, x -> vonNeumann(x, param.s^halfOfSystems, param.s^(param.numOfSys - halfOfSystems)))[1])
		end
		res
	end
	function makeXProjectors(s::Int64)
		out::Array{Array{Complex{Float64},2},1} = []
		a = lowOp(s)
		ad = copy(a')
		X = a + ad
		ev = eigvecs(X)
		for i in 1:s
			v = ev[:, i]
			normalize!(v)
			push!(out, v * v')
		end
		return out
	end
end

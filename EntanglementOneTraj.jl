using Plots, MIPTM, ParametersModule, Profile
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function firstState(sp)
	L = sp.L; N = sp.N; dim = sp.dim; cap = sp.cap
	basisState::Array{Int64,1} = first_state(L, N, cap)
	state = zeros(dim)
	state[find_index(basisState, cap)] = 1.
	return state
end

function halfL_BosonNumbermean(Ψ, L, N, dim; cap=N) # Numbers of bosons in the half of the chain
	nₕ = number(L, N, dim, 1; cap=cap)
	for i in 2:Int(round(L/2))
		nₕ .+= number(L, N, dim, i; cap=cap)
	end
	return expVal(Ψ, nₕ)
end
function halfL_BosonNumbermean(Ψ, sp)
	return halfL_BosonNumbermean(Ψ, sp.L, sp.N, sp.dim; cap=sp.cap)
end
function l_BosonNumbermean(Ψ, L, N, dim; cap=N) # Number of bosons at l in the middle of the chain
	nₕ = number(L, N, dim, round(Int64,L/2); cap=cap)
	#nₕ = number(L, L)
	return expVal(Ψ, nₕ)
end
function l_BosonNumbermean(Ψ, sp)
	return l_BosonNumbermean(Ψ, sp.L, sp.N, sp.dim; cap=sp.cap)
end

function f(prob, L)
	entForState(Ψ, p) = entanglement_entropy(p.sp.L, p.sp.N, Ψ, Int(p.sp.L / 2); cap=p.sp.cap)
	entanglement(sol, p) = calcMean(sol, Ψ -> entForState(Ψ, p))
	halfLNexpval(sol, p) = calcMean(sol, Ψ -> halfL_BosonNumbermean(Ψ, p.sp))
	lNexpval(sol, p) = calcMean(sol, Ψ -> l_BosonNumbermean(Ψ, p.sp))

	N = Int(L / 2);
	sp = SystemParameters(L=L, N=N)
	state = oneZeroState(sp)
	projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=prob, f=100.0, projOp=projOp)
	@time bhp = BoseHubbardParameters(sp=sp, U=5.0, J=1.0)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=6, dt=0.02, time=40.0, traj=1000)
	@time sol = MIPT(p, projectAfterTimeStep = true)
	println("ent")
	@time res = entanglement(sol, p)
	#@time res = calcMean(sol, Ψ->entanglement_entropy(p.sp.L, p.sp.N, Ψ, Int(L / 2), cap = p.sp.cap))
	#popfirst!(res)
	#popfirst!(p.t.times)
	plot(p.t.times, res, label="ent", title = "L = $(L), p = $(p.pp.p)", legend = :outertopright)
	res = halfLNexpval(sol, p)
	plot!(p.t.times, res, label="half chain n")
	res = lNexpval(sol, p)
	plot!(p.t.times, res, label="middle site n")
	savefig("L$(L)_p$(p.pp.p).png")
end
#BLAS.set_num_threads(1)
@time let 
	for i in [4,6,8]
		f(0.01, i)
		#f(0.1, i)
	end
println("end")
end
#f(0.06)
#f(0.16)

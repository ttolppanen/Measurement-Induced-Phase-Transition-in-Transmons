using Plots, MIPTM, ParametersModule
include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function firstState(sp)
	L = sp.L; N = sp.N; dim = sp.dim; cap = sp.cap
	basisState::Array{Int64,1} = first_state(L, N, cap)
	state = zeros(dim)
	state[find_index(basisState, cap)] = 1.
	return state
end

function f(prob)
	L = 6; N = L;
	sp = SystemParameters(L=L, N=N)
	state = firstState(sp)
	projOp = generateProjectionOperators(sp)
	pp = ProjectionParameters(p=prob, f=100.0, projOp=projOp)
	@time bhp = BoseHubbardParameters(sp=sp, U=5.0)
	p = Parameters(sp=sp, pp=pp, bhp=bhp, Ψ₀=state, sdim=3, dt=0.02, time=50.0, traj=100)
	@time sol = MIPT(p, projectAfterTimeStep = true)
	@time res = calcMean(sol, Ψ->entanglement_entropy(p.sp.L, p.sp.N, Ψ, 3, cap = p.sp.cap))
	#popfirst!(res)
	#popfirst!(p.t.times)
	plot(p.t.times, res, label="p=$(p.pp.p)", legend=:left)
end

@time f(0.0)
#f(0.06)
#f(0.16)

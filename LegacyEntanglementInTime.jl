using MIPTMLegacy
using Plots

function f()
	Ψ₀ = [[0.0, 0+0im, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]
	measOp = makeXProjectors(3)
	p = ParametersConstructor(Ψ₀=Ψ₀, measOp=measOp, t=(0.0, 0.02, 40.0),
	Γ=0.16, ω=0.0, U=0.14, J=1.0, traj=100)
	sol = MIPT(p)
	#res = calcMean(sol, x->vonNeumannHalfOfSystem(x, p))
	res2 = calcMean(sol, x->expVal(x, p.op.nAll))
	#plot!(p.t.times, res, label=p.Γ)
	plot!(p.t.times, res2, label="")
end

@time f()

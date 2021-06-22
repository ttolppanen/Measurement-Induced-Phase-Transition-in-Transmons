using MIPTM
using Plots

function f()
	Ψ₀ = [[1.0, 0im], [0, 1], [1, 0], [0, 1], [1, 0], [0, 1]]
	measOp = [[1.0 0im; 0 0], [0 0 ; 0 1]]
	p = ParametersConstructor(Ψ₀=Ψ₀, measOp=measOp, t=(0.0, 0.1, 10.0),
	Γ=0.2, ω=0.0, U=1.0, J=1.0, traj=10)
	sol = MIPT(p)
	res = calcMean(sol, x -> vonNeumann(x, p.s^2, p.s^2))
	res1 = calcMean(sol, x->expVal(x, p.op.nAll))
	plot(p.t.times, res, ylim=[0,2])
	plot!(p.t.times, res1)
end

@time f()

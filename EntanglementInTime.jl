using MIPTM
using Plots

function f()
	Ψ₀ = [[1.0, 0im, 0], [0, 0, 1], [1, 0, 0], [1, 0, 0]]
	measOp = [[1.0 0im 0; 0 0 0; 0 0 0], [0 0im 0; 0 1 0; 0 0 0], [0 0im 0; 0 0 0; 0 0 1]]
	p = ParametersConstructor(Ψ₀=Ψ₀, measOp=measOp, t=(0.0, 0.01, 10.0),
	Γ=0.2, ω=0.0, U=10.0, J=1.0, traj=1)
	sol = MIPT(p)
	res = calcMean(sol, x->vonNeumannHalfOfSystem(x, p))
	res2 = calcMean(sol, x->expVal(x, p.op.n𝐼[2]))
	res3 = calcMean(sol, x->expVal(x, p.op.n𝐼[3]))
	plot(p.t.times, res, ylim=[0,2])
	plot!(p.t.times, res2)
	plot!(p.t.times, res3)
end

@time f()

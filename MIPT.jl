using MIPTM
using Plots

function f()
	Ψ₀ = [[0.0+0im, 1.0, 0.0, 0.0], [0.0+0im, 1.0, 0.0, 0.0], [0.0+0im, 0.0, 1.0, 0.0]]
	measOp = [[1 0 0 0; 0 0 0 0.0im; 0 0 0 0; 0 0 0 0], [0 0 0 0; 0 1 0 0.0im; 0 0 0 0; 0 0 0 0],
	[0 0 0 0; 0 0 0 0.0im; 0 0 1 0; 0 0 0 0], [0 0 0 0; 0 0 0 0.0im; 0 0 0 0; 0 0 0 1]]
	p = Parameters(Ψ₀=Ψ₀, measOp=measOp, t=(0.0, 0.01, 5.0), p=0.5, f=1.0, ω=1.0, U=1.0, J=1.0, traj=100)
	sol = MIPT(p)
	res1 = calcMean(sol, x->expVal(x, p.op.n𝐼[1]))
	plot(p.t.times, res1, ylims=(0,3))
end

f()

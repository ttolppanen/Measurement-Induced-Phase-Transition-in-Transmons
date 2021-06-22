using MIPTM
using DifferentialEquations
using Plots
using Random
using LinearAlgebra

function f()
	Ψ₀ = [[1.0+0im, 0.0, 0.0, 0.0], [0.0+0im, 1.0, 0.0, 0.0], [0.0+0im, 0.0, 1.0, 0.0]]
	measOp = [[1 0 0 0; 0 0 0 0.0im; 0 0 0 0; 0 0 0 0], [0 0 0 0; 0 1 0 0.0im; 0 0 0 0; 0 0 0 0],
	[0 0 0 0; 0 0 0 0.0im; 0 0 1 0; 0 0 0 0], [0 0 0 0; 0 0 0 0.0im; 0 0 0 0; 0 0 0 1]]
	p = Parameters(Ψ₀=Ψ₀, measOp=measOp, t=(0.0, 0.01, 5.0), Γ=0.5, ω=1.0, U=1.0, J=1.0)
	sol = MIPT(p)
	res1 = calcMean(sol, x->expVal(x, p.op.n𝐼[1]))
	res2 = calcMean(sol, x->expVal(x, p.op.n𝐼[2]))
	res3 = calcMean(sol, x->expVal(x, p.op.n𝐼[3]))
	res4 = calcMean(sol, x->expVal(x, p.op.nAll))
	plot(p.t.times, res1)
	plot!(p.t.times, res2)
	plot!(p.t.times, res3)
	plot!(p.t.times, res4)
end

f()

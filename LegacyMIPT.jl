using MIPTMLegacy, LinearAlgebra
using Plots

function f()
	Ψ₀ = [[0.0, 1 + 0im, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0]]
	measOp = makeXProjectors(3)
	rates = 0.02:0.02:0.18
	p = ParametersConstructor(Ψ₀=Ψ₀, measOp=measOp, t=(0.0, 0.02, 10.0),
	Γ=0.01, ω=0.0, U=4.0, J=1.0, traj=100)
	res = entanglementAndMeasProbability(p, rates)
	plot!(rates, res, label="L = 6")
end

f()

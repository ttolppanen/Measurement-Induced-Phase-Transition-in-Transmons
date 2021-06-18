using MIPTM
using Plots

function f()
	Ψ₀ = [[0.0+0im, 1.0, 0.0], [0.0+0im, 1.0, 0.0], [0.0+0im, 1.0, 0.0]]
	measOp = [[1 0.0 0; 0im 0 0; 0 0 0], [0 0 0; 0 1 0; 0 0 0], [0 0 0; 0 0 0; 0 0 1]]
	probabilities = 0.0:0.1:0.4
	p = ParametersConstructor(Ψ₀=Ψ₀, measOp=measOp, t=(0.0, 0.5, 20.0), p=0.0, f=10.0, ω=1.0, U=1.0, J=1.0, traj=10)
	res = entanglementAndMeasProbability(p, probabilities)
	plot!(probabilities, res)
end

@time f()

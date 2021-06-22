using MIPTM
using Plots

function f()
	Ψ₀ = [[1.0, 0im], [0, 1], [1, 0], [0, 1]]
	measOp = [[1.0 0im; 0 0], [0 0 ; 0 1]]
	probabilities = 0.05:0.05:0.2
	p = ParametersConstructor(Ψ₀=Ψ₀, measOp=measOp, t=(0.0, 0.1, 10.0),
	p=0.0, f=1.0, ω=1.0, U=1.0, J=1.0, traj=100)
	res = entanglementAndMeasProbability(p, probabilities)
	plot(probabilities, res, ylim=[0, 1])
end

@time f()

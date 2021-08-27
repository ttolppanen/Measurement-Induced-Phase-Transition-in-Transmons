using Plots

function f()
	prob = 0.0:0.01:0.09
	plot(prob, x->1.0/(x^2 + 3))
	plot!(prob, x->0)
end

f()

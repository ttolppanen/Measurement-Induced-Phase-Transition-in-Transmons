include.(["OllisCode/Operators.jl", "OllisCode/Time.jl", "OllisCode/Density.jl",
		"OllisCode/Basis.jl", "OllisCode/Entropy.jl"])

function f()
	L = 3
	N = 3
	p1 = projector(L, N, 1, 1)
	display(p1)
	singleSiteOperator = [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0]
	p2 = generalizeSingleSiteOperator(L, N, 1, singleSiteOperator)
	display(p2)
end
f()

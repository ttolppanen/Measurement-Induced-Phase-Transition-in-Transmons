using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Basis.jl"])

function f()
	print(dimension(3, 10))
end

f()

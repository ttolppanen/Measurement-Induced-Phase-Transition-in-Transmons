using MIPTM, Plots
include("OllisCode/Basis.jl")


function f()
	L = 8
	N = L
	cap = N
	d = dimensions(L, N, cap)
	display(d)
end

f()

using MIPTM, Plots
include("OllisCode/Basis.jl")


function f()
	L = 8
	N = L
	cap = 4
	d = dimensions(L, N, cap=cap)
	display(d)
end

f()

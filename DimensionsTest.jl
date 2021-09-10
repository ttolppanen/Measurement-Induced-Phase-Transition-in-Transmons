using MIPTM, Plots
include("OllisCode/Basis.jl")


function f()
	L = 20
	N = L
	cap = 20
	display("ASD")
	@time d = dimensions(L, N, cap=cap)
	display(d)
end

f()

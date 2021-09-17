using MIPTM, Plots
include("OllisCode/Basis.jl")

function calcD(L, N, cap)
	cut = find_index(first_state(L, N, cap), N) - 1
	d = dimensions(L, N)
	return d - cut
end

function f()
	L = 4
	N = L
	cap = 3
	@time dReal = dimensions(L, N, cap=cap)
	dMaybe = calcD(L, N, cap)
	display(dReal)
	display(dMaybe)
end
function a()
	L = 4
	N = L
	s = first_state(L, N, N)
	display("HOOI")
	display(s)
	d = dimensions(L, N)
	for _ in 2:d
		next!(s, N)
		display(s)
	end
end

a()

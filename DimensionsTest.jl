using MIPTM, Plots
include("OllisCode/Basis.jl")

function calcD(L, N, cap)
	cut = find_index(first_state(L, N, cap), N) - 1
	d = dimensions(L, N)
	return d - cut
end
function findUpperN(L, d)
	N = 0
	L = Int((L + 2) / 2)
	while true
		N += 1
		dim = dimensions(L, N)
		if dim == d
			display(N)
			break
		elseif dim > d
			display("EII")
			break
		end
	end
end
function totalDim(L, N)
	return N^L
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
	L = 8
	N = 5
	d = dimensions(L, N)
	display(d)
	display(totalDim(4, 5))
end

a()

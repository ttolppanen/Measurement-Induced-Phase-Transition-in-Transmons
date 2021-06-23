using Plots, MIPTM
include.(["OllisCode/Operators.jl", "OllisCode/Basis.jl"])

function makeN1(L, N)
	d = dimension(L, N)
	n = spzeros(d, d)
	basis_vector = zeros(Int64, L)
	basis_vector[1] = N
	for i in 1:d
		n[i, i] = basis_vector[1]
		next!(basis_vector)
	end
	n
end

function f()
	n = makeN1(3, 2)
	display(n)
end

f()

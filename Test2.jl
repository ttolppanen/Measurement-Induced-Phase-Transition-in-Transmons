using LinearAlgebra

function expM(M)
	va, U = eigen(M)
	d = U' * M * U
	d .= exp(Diagonal(d))
	return U * d * U'
end

f()

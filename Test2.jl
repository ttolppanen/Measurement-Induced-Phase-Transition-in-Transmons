using LinearAlgebra, Arpack

function fs(s, i)
	return i * s
end

function f()
	m = spzeros(3, 3)
	m[1,1] = 2
	m[2,2] = 3
	m[3,3] = 4
	Z = svds(m)
	display(exp.(m))
	display()
end

f()

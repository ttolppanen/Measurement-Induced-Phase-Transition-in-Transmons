using MIPTM, Plots

function f()
	v = [1,0,1,0,1]
	L = length(v)
	N = sum(v)
	d = binomial(L, N)
	displayVec(v)
	for i in 1:d
		q_next!(v)
		displayVec(v)
	end
end

function displayVec(v)
	i = q_find_index(v)
	display(i)
	show(v)
end

f()

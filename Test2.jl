using Plots

function asdasd(mat, dim)
	for i in 1:dim
		for j in 1:dim
			mat[j, i] = 0
		end
	end
end

function f()
	w = [1, 2, 3]
	m1 = [1 2 3; 3 4 5; 1 2 3]
	m2 = [1,2,3]
	@time a = w' * w
end

f()

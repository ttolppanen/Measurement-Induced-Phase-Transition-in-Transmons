using Plots

function asdasd(mat, dim)
	for i in 1:dim
		for j in 1:dim
			mat[j, i] = 0
		end
	end
end

function f()
	dim = 10000
	a = ones(dim, dim)
	a .= zeros(dim, dim)
end

@time f()

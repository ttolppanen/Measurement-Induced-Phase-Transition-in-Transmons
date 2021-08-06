using Combinatorics

function allPartitions(N, cap)
	out = []
	p = collect(partitions(N))
	for pᵢ in p
		if length(pᵢ[pᵢ .> cap]) == 0
			tempOut = []
			for i in 1:cap
				push!(tempOut, length(pᵢ[pᵢ .== i]))
			end
			push!(out, tempOut)
		end
	end
	return out
end
function dimensions(L, N, cap)
	d = 0
	p = allPartitions(N, cap)
	for pᵢ in p
		dAdd = factorial(L) / factorial(L - sum(pᵢ))
		for i in pᵢ
			dAdd /= factorial(i)
		end
		d += dAdd
	end
	return Int(d)
end

function f()
	display(dimensions(3, 2, 2))
end

@time f()

using Plots
using BSON: @load

function readBsonFile(title)
	path = pwd() * "/Plots/" * title * "/data.bson"
	@load path x y
	prob = x
	mean = y[1]
	stantardDeviation =  [sqrt.(i) for i in y[2]]
	return prob, mean, stantardDeviation
end

function f()
	p = 0.06
	v = 5.0
	C = 2.00
	prob, ent = readBsonFile("NTest")
	x = []
	y = []
	L = [21, 56, 126, 252, 462, 792, 1287]
	for i in 1:length(L)
		push!(x, [(j - p)*L[i]^(1/v) for j in prob])
		push!(y, L[i]^(-C/v) .* ent[i])
	end
	plot(x, y)
end

f()

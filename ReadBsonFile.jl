using BSON: @load

function readBsonFile(title)
	path = pwd() * "/Plots/" * title * "/data.bson"
	@load path x y
	prob = x
	mean = y[1]
	stantardDeviation =  [sqrt.(i) for i in y[2]]
	return prob, mean, stantardDeviation
end

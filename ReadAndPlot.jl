using Plots
using BSON: @load

function readAndPlot(title)
	path = pwd() * "/Plots/" * title * "/data.bson"
	@load path x y
	display(typeof(y))
	plot(x, y)
end
function readAndPlot!(title)
	path = pwd() * "/Plots/" * title * "/data.bson"
	@load path x y
	plot!(x, y)
end

function f()
	readAndPlot("ELV_S_5000_1000_300")
end

f()

using Plots
using BSON: @load

function readAndPlot(title)
	path = pwd() * "/Plots/" * title * "/data.bson"
	@load path x y
	plot(x, y, color="blue")
end
function readAndPlot!(title)
	path = pwd() * "/Plots/" * title * "/data.bson"
	@load path x y
	plot!(x, y, color="red")
end

function f()
	readAndPlot("ELV_N_L_05")
	readAndPlot!("DisorderTest")
end

f()

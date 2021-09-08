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
	readAndPlot("Fixed_ELV_S_d3_dis_Attractive_Insulator_10000_5000_300")
	#readAndPlot!("ELV_S_d3_dis_Attractive_Insulator_7000_2000_100")
	readAndPlot!("ELV_S_2020test_5000_2000_100")
end

f()

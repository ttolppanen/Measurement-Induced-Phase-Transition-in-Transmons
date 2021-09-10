import fssa
import numpy as np
import matplotlib.pyplot as plt
import julia
j = julia.Julia()
j.include("ReadBsonFile.jl")
from julia import Main


def findBestValue(p, S, dS):
	bestQuality = 1000
	bestRet = 0
	loopValues = np.arange(0.5, 2.5, 0.3)
	for i in loopValues:
		for j in loopValues:
			ret = fssa.autoscale([4, 6, 8], p, S, dS, 0.07, i, j)
			sd = fssa.scaledata([4, 6, 8], p, S, dS, ret.rho, ret.nu, ret.zeta)
			quality = fssa.quality(sd.x, sd.y, sd.dy)
			if quality < bestQuality:
				bestRet = ret
				bestQuality = quality
	return bestRet


p, S, dS = Main.readBsonFile("Fixed_ELV_S_d3_dis_Attractive_Insulator_10000_5000_300")
# ret = fssa.autoscale([4, 6, 8], p, S, dS, 0.06, 1.6, 0.5)
ret = findBestValue(p, S, dS)
sd = fssa.scaledata([4, 6, 8], p, S, dS, ret.rho, ret.nu, ret.zeta)
plt.plot(sd.x.T, sd.y.T)
plotName = "Insulator"
plotName += "Min0"
plotName += "_p{:.3f}_dp{:.3f}".format(ret.rho, ret.drho)
# plt.savefig(plotName + ".png")
print(ret)
print("quality = " + str(fssa.quality(sd.x, sd.y, sd.dy)))
plt.show()

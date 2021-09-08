import fssa
import numpy as np
import matplotlib.pyplot as plt
import julia
j = julia.Julia()
j.include("ReadBsonFile.jl")
from julia import Main


def findBestValue(p, S, dS):
	leastError = 1000
	bestRet = 0
	loopValues = np.arange(0.1, 2.0, 0.1)
	for i in loopValues:
		for j in loopValues:
			ret = fssa.autoscale([4, 6, 8], p, S, dS, 0.06, i, j)
			error = sum(ret.errors)
			if error < leastError:
				bestRet = ret
				leastError = error
	return bestRet


p, S, dS = Main.readBsonFile("Fixed_ELV_S_d3_dis_Attractive_Superfluid_7000_2000_100")
# ret = fssa.autoscale([4, 6, 8], p, S, dS, 0.06, 2.0, 0.5)
ret = findBestValue(p, S, dS)
auto_scaled_data = fssa.scaledata([4, 6, 8], p, S, dS, ret.rho, ret.nu, ret.zeta)
plt.plot(auto_scaled_data.x.T, auto_scaled_data.y.T)
print(ret)
plt.show()

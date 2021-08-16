using MIPTM, Plots


function f()
	L = 10; N = Int(L/2);
	cap = 6
	state = zeroOneState(L, N, cap)
	measOp = generateProjectionOperators(L, N, cap)
	p = ParametersConstructor(L=L, N=N, cap=cap, sdim=3, measOp=measOp, traj=1, dt=0.02, time=3.0, p=0.1, f=1.0, U=0.14, J=1.0, Ψ₀=state)
	prop = propagator(p.𝐻, p.Ψ₀, p.sdim, p.t.dt)
	pic = plot(abs2.(prop))
	display(pic)
end

f()

using MIPTM
using DifferentialEquations
using Plots

function schrodinger(Ψ, p, t)
	-1im * p.𝐻 * Ψ
end
function f()
	p = Parameters(numOfSys=3, s=3, t=(0.0, 0.01, 5.0), p=0.5, f=1.0, ω=1.0, U=5.0, J=1.0)
	Ψ₀ = kronForMany([[0.0+0im, 1.0, 0.0], [0.0+0im, 1.0, 0.0], [0.0+0im, 0.0, 1.0]])
	prob = ODEProblem(schrodinger, Ψ₀, p.t.Δt, p, saveat = p.t.dt)
	enProb = EnsembleProblem(prob, safetycopy=true)
	sol = solve(enProb, trajectories = p.traj)
	sol = ensSolToList(sol, p.t.times)
	res1 = calcMean(sol, x->expVal(x, p.op.n𝐼[1]))
	res2 = calcMean(sol, x->expVal(x, p.op.n𝐼[2]))
	res3 = calcMean(sol, x->expVal(x, p.op.n𝐼[3]))
	res4 = calcMean(sol, x->expVal(x, p.op.nAll))
	plot(p.t.times, res1)
	plot!(p.t.times, res2)
	plot!(p.t.times, res3)
	plot!(p.t.times, res4)
end

f()

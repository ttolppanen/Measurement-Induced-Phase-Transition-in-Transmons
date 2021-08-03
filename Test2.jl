using MIPTM

function ELVSavefig(;title::String, p::Parameters, initialState::String)
	folderName = "ELV_" * title
	path = pwd() * "/Plots/Last Values of Entanglement/" * folderName
	mkpath(path)
	io = open(path * "/data.txt", "w")
	println(io, "L = " * string(p.L))
	println(io, "N = " * string(p.N))
	println(io, "sdim = " * string(p.sdim))
	println(io, "p = " * string(p.p))
	println(io, "f = " * string(p.f))
	println(io, "dt = " * string(p.t.dt))
	println(io, "t = " * string(p.t.times[end]))
	println(io, "Ψ₀ = " * initialState)
	close(io)
	savefig(pwd() * "/Plots/Last Values of Entanglement/" * folderName * "/" * folderName * ".png")
end
function f()
	L = 3
	N = L
	state = zeroOneState(L, N)
	measOp = generateProjectionOperators(L, N)
	p = ParametersConstructor(L=L, N=N, sdim=10, measOp=measOp, traj=100, dt=0.02, time=10.0, p=0.0, f=1.0, U=0.14, J=1.0, Ψ₀=state)
	ELVSavefig(p=p, title="C_1000_100_10", initialState="101010")
end

f()

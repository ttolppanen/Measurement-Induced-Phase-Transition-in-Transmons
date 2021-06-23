using Plots, MIPTM
include.(["Operators.jl", "Time.jl", "Density.jl", "Basis.jl"])

function makeN1(L, N)
	d = dimension(L, N)
	n = spzeros(d, d)
	basis_vector = zeros(Int64, L)
	basis_vector[1] = N
	for i in 1:d
		n[i, i] = basis_vector[1]
		next!(basis_vector)
	end
	n
end

function f()
	L = 8; N = 4;
	HU, HJ = split_hamiltonian(L, N)
	H = 0.5 .* HU .+ HJ
	basis_vector = zeros(Int64, L)
	basis_vector[1] = N
	state = zeros(dimension(L, N))
	state[find_index(basis_vector)] = 1.;

	dt = 0.1; sdim = 10; steps = 100;
	out = []
	n = makeN1(L, N)
	push!(out, real(state' * n * state))
	for i in 1:steps
		state = propagate(H, state, sdim, dt)
		push!(out, real(state' * n * state))
	end

	time = collect(range(0, 10, length = steps + 1))
	plot(time, out)
end

@time f()

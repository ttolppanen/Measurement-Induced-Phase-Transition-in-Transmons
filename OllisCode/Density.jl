using LinearAlgebra
include("Basis.jl")

function density(L, N, state)
    density = zeros(N + 1, L)
    basis_vector = zeros(Int64, L)
    basis_vector[1] = N
    dim = dimension(L, N)
    for i in 1:dim
        if i != 1
            next!(basis_vector)
        end

        #index = find_index(basis_vector)
        for j in 1:L
            density[basis_vector[j] + 1, j] += abs2(state[i])
        end
    end

    return density
end


function density_profile(L, N, state)
    density = zeros(L)
    basis_vector = zeros(Int64, L)
    basis_vector[1] = N
    dim = dimensionOlli(L, N)
    for i in 1:dim
        if i != 1
            nextOlli!(basis_vector)
        end

        #index = find_index(basis_vector)
        for j in 1:L
            density[j] += abs2(state[i]) * basis_vector[j]
        end
    end

    return density
end


function print_density_profile(L, N, rho)
    rhop = zeros(L)
    ns = 1. .* collect(1:N)
    for j in 1:L
        rhop[j] = sum(rho[2:end, j]' * ns)
    end

    display(rhop')
    println()
end


function local_energy(rho, dis = 0)
    E = 0.
    N, L = size(rho)
    N -= 1
    Ns = 1. .* collect(0:N)
    NNs = Ns .* (Ns .- 1.)
    for j in 1:L
        if dis != 0
            E += dis[j] * sum(rho[:, j] .* Ns)
        end

        E -= 0.5 * sum(rho[:, j] .* NNs)
    end

    return E
end

using LinearAlgebra, SparseArrays, KrylovKit
include.(["Basis.jl", "Operators.jl"])

function entanglement_entropy(L, N, state, sites)
    sub_dim = dimension_open(sites, N)
    schmidt = zeros(ComplexF64, sub_dim, dimension_open(L - sites, N))
    basis_vector = zeros(Int64, L)
    basis_vector[1] = N
    for i in 1:dimensionOlli(L, N)
        if i != 1
            nextOlli!(basis_vector)
        end

        schmidt[find_index_open(basis_vector[1:sites]), find_index_open(basis_vector[sites + 1: end])] = state[i]
    end

    S = svd(schmidt).S
    entropy = 0.0
    for i in 1:sub_dim
        if S[i] > 1e-16
            alpha2 = abs2(S[i])
            entropy += -alpha2 * log(alpha2)
        else
            break
        end
    end

    return entropy
end



function test()
    L = 8; N = 4
    HU, HJ = split_hamiltonian(L, N)
    H = 0.05 .* HJ .+ HU
    vals, vecs, info = eigsolve(H, 1, :SR)
    state = vecs[1]
    println(entanglement_entropy(L, N, state, 4))
    println("should be bit above ln 2 = 0.693...")
end

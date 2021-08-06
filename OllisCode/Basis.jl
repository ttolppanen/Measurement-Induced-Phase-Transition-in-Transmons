using Combinatorics

function isThereTooManyBosons(v, cap::Int64)
    for i in v
        if i > cap
            return true
        end
    end
    return false
end
function next!(v, cap)
    N = sum(v)
    nextOlli!(v)
    if N != cap && isThereTooManyBosons(v, cap)
        next!(v)
    end
end
function find_index(v, cap)
    N = sum(v)
    L = length(v)
    if N == cap
        return find_indexOlli(v)
    else
        loopState = first_state(L, N, cap)
        d = binomial(L, N)
        for i in 1:d
            if v == loopState
                return i
            end
            q_next!(loopState)
        end
    end
end
function first_state(L, N, cap)::Vector{Int} #cap = max amount of bosons on one site...
    fullStates = Int(floor(N/cap))
    bosonsLeftOver = N % cap
    return [cap * ones(fullStates); bosonsLeftOver; zeros(L - fullStates - 1)]
end
function allPartitions(N, cap)
    out = []
    p = collect(partitions(N))
    for pᵢ in p
        if length(pᵢ[pᵢ .> cap]) == 0
            tempOut = []
            for i in 1:cap
                push!(tempOut, length(pᵢ[pᵢ .== i]))
            end
            push!(out, tempOut)
        end
    end
    return out
end
function dimensions(L, N, cap)
    if cap == N
        return dimensionOlli(L, N)
    else
        d = 0
        p = allPartitions(N, cap)
        for pᵢ in p
            dAdd = factorial(L) / factorial(L - sum(pᵢ))
            for i in pᵢ
                dAdd /= factorial(i)
            end
            d += dAdd
        end
        return Int(d)
    end
end

function nextOlli!(vector)
    L = length(vector)
    N = sum(vector)
    nk = 1
    for j in reverse(1:L - 1)
        if vector[j] != 0
            nk = j
            break
        end
    end

    if vector[nk] > 0
        vector[nk] -= 1
        vector[nk + 1] = N - sum(vector[1:nk])
        vector[nk + 2:end] *= 0
    else
        vector[:] *= 0
        vector[1] = N
    end
end


function dimensionOlli(L, N)::Int
    d = 1.
    for i in 1:min(N, L - 1)
        d *= (L + N - i) / i
    end

    return round(Int64, d)
end


function dimension_open(L, N)::Int
    dim = 1
    for n in 1:N
        dim += dimensionOlli(L, n)
    end

    return dim
end


function find_indexOlli(vector::Vector{Int})::Int
    L = size(vector, 1)
    N = sum(vector)
    index = 1

    for k in 1:L - 1
        index += binomial(N + L - 1 - k - sum(vector[1:k]), L - k)
    end

    return index
end


function find_index_open(vector::Vector{Int})::Int64
    L = size(vector, 1)
    N = sum(vector)
    if N > 0
        index = 2
        for i in 1:N - 1
            index += dimensionOlli(L, i)
        end

        for k in 1:L - 1
            index += binomial(N + L - 1 - k - sum(vector[1:k]), L - k)
        end

        return index
    else

        return 1
    end
end


function print_state(L, N, state, cutoff = 0.99)
    dim = dimension(L, N)
    basis_vector = zeros(Int64, L)
    basis_vector[1] = N
    sorting_thing = []
    basis = []
    for i in 1:dim
        push!(sorting_thing, abs(state[i]))
        fock = "|" * prod(string.(basis_vector)) * "⟩"
        push!(basis, fock)
        next!(basis_vector)
    end

    indices = sortperm(sorting_thing, rev = true)
    total = 0.
    for i in 1:dim
        total += abs2(state[indices[i]])
        println(round(state[indices[i]], digits = 4), " ", basis[indices[i]], "         ", round(total, digits = 4))
        if total > cutoff
            break
        end
    end

    println()
end

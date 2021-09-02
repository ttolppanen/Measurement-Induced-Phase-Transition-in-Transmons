using LinearAlgebra, SparseArrays, ParametersModule

function lanczos(H, state, sdim)
    dim = length(state)
    h = zeros(ComplexF64, sdim, sdim)
    V = zeros(ComplexF64, dim, sdim)
    V[:, 1] .= normalize(state)
    w = H * state
    h[1, 1] = w' * state
    w .-= h[1, 1] * state
    for j in 2:sdim
        beta = norm(w)
        V[:, j] .= w / beta
        w = H * V[:, j]
        h[j, j] = w' * V[:, j]
        w .-= (h[j, j] * V[:, j] + beta * V[:, j - 1])
        h[j - 1, j] = beta
        h[j, j - 1] = beta
    end

    return V, h
end
function lanczos!(H, state, dim, sdim, V, h, w)
    V .= zeros(ComplexF64, dim, sdim)
    h .= zeros(ComplexF64, sdim, sdim)
    V[:, 1] .= normalize(state)
    w .= H * state
    h[1, 1] = w' * state
    w .-= h[1, 1] * state
    for j in 2:sdim
        beta = norm(w)
        V[:, j] .= w / beta
        w = H * V[:, j]
        h[j, j] = w' * V[:, j]
        w .-= (h[j, j] * V[:, j] + beta * V[:, j - 1])
        h[j - 1, j] = beta
        h[j, j - 1] = beta
    end
end


function arnoldi(H, state, sdim)
    dim = length(state)
    h = zeros(ComplexF64, sdim, sdim)
    V = zeros(ComplexF64, dim, sdim)
    v = normalize(state)
    V[:, 1] .= v
    for j in 1:sdim
        v .= H * v
        for i in 1:j
            h[i, j] = V[:, i]' * v
            v .-= h[i, j] .* V[:, i]
        end

        if j != sdim
            h[j + 1, j] = norm(v)
            v ./= h[j + 1, j]
            V[:, j + 1] .= v
        end
    end

    return V, h
end


function propagate!(p::Parameters, H::SparseMatrixCSC{Float64,Int64}, state)
    propagate!(H, state, p.t.dt, p.sp.dim, p.sdim, p.V[Threads.threadid()], p.h[Threads.threadid()], p.w[Threads.threadid()])
end
function propagate!(H::SparseMatrixCSC{Float64,Int64}, state, dt::Float64, dim::Int64, sdim::Int64, V, h, w)
    lanczos!(H, state, dim, sdim, V, h, w)
    state .= normalize(V * exp(-im * dt * h)[1, :])
end
function propagate!(H::SparseMatrixCSC{Float64,Int64}, state, sdim::Int64, dt; algo = lanczos, normalise = true)
    V, h = algo(H, state, sdim)
    if normalise == true
        state .= normalize(V * exp(-im * dt * h)[1, :])
    else
        state .= V * (exp(-im * dt * h)[1, :])
    end
end


function propagator(H, state, sdim, dt)
    #The effective Hamiltonian is time-independent, and the exponent needs to be computed only once
    #NB this doesn't work for imaginary time evolution, also seems to cause a small cumulative error
    dim = length(state)
    h = zeros(ComplexF64, sdim, sdim)
    V = zeros(ComplexF64, dim, sdim)
    V[:, 1] .= normalize(state)
    w = H * state
    h[1, 1] = w' * state
    w .-= h[1, 1] * state
    for j in 2:sdim
        beta = norm(w)
        V[:, j] .= w / beta
        w = H * V[:, j]
        h[j, j] = w' * V[:, j]
        w .-= (h[j, j] * V[:, j] + beta * V[:, j - 1])
        h[j - 1, j] = beta
        h[j, j - 1] = beta
    end


    return exp(-im * dt * h)[1, :]
end


function propagate(H, state, propagator::Vector)
    #If the effective Hamiltonian has been computed beforehand, only three Krylov vectors need to be saved at a time.
    krylov1 = copy(state)
    krylov2 = H * krylov1
    krylov2 .-= (krylov1' * krylov2) .* krylov1
    normalize!(krylov2)
    state .*= propagator[1]
    state .+= propagator[2] .* krylov2
    for i in 3:length(propagator)
        krylov3 = H * krylov2
        krylov3 .-= (krylov1' * krylov3) .* krylov1
        krylov3 .-= (krylov2' * krylov3) .* krylov2
        normalize!(krylov3)
        state .+= propagator[i] .* krylov3
        if i < length(propagator)
            krylov1 = copy(krylov2)
            krylov2 = copy(krylov3)
        end
    end


    return normalize(state)

end




#~ function toeplitz(vals)
#~     dim = Integer(floor((length(vals) + 1) / 2))
#~     T = zeros(ComplexF64, dim, dim)
#~     vals = reverse(vals)
#~     for i in 1:dim
#~         T[i, :] = vals[end - i - dim + 1:end - i]
#~     end

#~     return T
#~ end


#~ using Polynomials, ToeplitzMatrices


#~ function moments(H, state, n)
#~     dim = length(state)
#~     V = zeros(n, dim)
#~     V[1, :] = H * state
#~     for i in 2:n
#~         V[i, :] .= H * V[i - 1, :]
#~     end

#~     return V
#~ end


#~ function loschmidt_zeros(ms, nzs, time, eps = 1e-2)
#~     cs = copy(ms)
#~     ncs = length(ms)
#~     for n in 2:ncs
#~         for m in 1:n - 1
#~             cs[n] -= binomial(n - 1, m - 1) * cs[m] * ms[n - m]
#~         end
#~     end

#~     for n in 1:ncs
#~         cs[n] *= (-1)^(n - 1) / factorial(big(n - 1))
#~     end


#~     ks = cs[end - nzs + 1:end]
#~#~     T = toeplitz(cs[collect(ncs - 2 * nzs + 1:ncs) .- 1])
#~     T = Toeplitz(cs[collect(ncs - nzs:ncs - 1)], cs[collect(reverse(ncs - 2 * nzs + 1:ncs - nzs))])
#~     as = T \ ks


#~     lambdas = roots(Polynomial(push!(reverse(as), -1.)))
#~     V = zeros(ComplexF64, nzs, nzs)
#~     for i in 1:nzs
#~         V[i, :] .= lambdas.^(i - 1)
#~     end

#~     ds = V \ ks ./ lambdas .^ (ncs - nzs + 1)


#~     out::Array{ComplexF64, 1} = []
#~     for i in 1:nzs
#~         if abs(1 - real(ds[i])) < eps
#~             push!(out, time - 1im / lambdas[i])
#~         end
#~     end

#~     return out
#~ end

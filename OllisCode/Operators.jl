using SparseArrays, ParametersModule
include("Basis.jl")
#=
function generalizeSingleSiteOperator(L, N, l, singleSiteOperator) #l=which site
    d = dimension(L, N)
    totalSystemOperator = spzeros(d, d) #create a d x d sparse matrix
    ket = zeros(Int64, L) #|i><j|
    ket[1] = N
    for i in 1:d
        bra = zeros(Int64, L) #|i><j|
        bra[1] = N
        for j in 1:d
            si = ket[l] + 1; sj = bra[l] + 1;
            totalSystemOperator[i, j] = singleSiteOperator[si, sj]
            next!(bra)
        end
        next!(ket)
    end
    return totalSystemOperator
end
=#

function projector(L, N, dim, l, n; cap=N)
    P = spzeros(dim, dim)
    basis_vector = first_state(L, N, cap)
    for i in 1:dim
        if i != 1
           next!(basis_vector, cap)
        end

        if basis_vector[l] == n
            P[i, i] = 1.
        end
    end

    return P
end

function number(L, N, dim, site; cap=N)
    n = spzeros(dim, dim)
    basis_vector = first_state(L, N, cap)
    for i in 1:dim
        if i != 1
           next!(basis_vector, cap)
        end

        n[i, i] += basis_vector[site]
    end

    return n
end


function numbers(L, N)
    n_l::Array{SparseMatrixCSC{Float64, Int64}, 1} = []
    for l in 1:L
        push!(n_l, number(L, N, l))
    end

    return n_l
end


function numbers_squared(L, N)
    nn_l::Array{SparseMatrixCSC{Float64, Int64}, 1} = []
    for l in 1:L
        n = number(L, N, l)
        push!(nn_l, n * n)
    end

    return nn_l
end


function current(L, N, site)
    dim = dimension(L, N)
    curr = spzeros(ComplexF64, dim, dim)
    basis_vector = zeros(Int64, L)
    basis_vector[1] = N
    for i in 1:dim
        if i != 1
           next!(basis_vector)
        end

        modified_vector = copy(basis_vector)
        if basis_vector[site] > 0
            if site < L
                modified_vector[site] -= 1
                modified_vector[site + 1] += 1
                index = find_index(modified_vector)
                curr[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1])
                curr[index, i] = -curr[i, index]
            else
                modified_vector[L] -= 1
                modified_vector[1] += 1
                index = find_index(modified_vector)
                curr[i, index] = sqrt(basis_vector[L] * modified_vector[1])
                curr[index, i] = -curr[i, index]
            end
        end
    end

    return curr .* -im
end


function split_hamiltonian(L, N; periodic = false)
    dim = dimensionOlli(L, N)
    #HD = spzeros(dim, dim)
    HJ = spzeros(dim, dim)
    HU = spzeros(dim, dim)
    basis_vector = zeros(Int64, L)
    basis_vector[1] = N
    for i in 1:dim
        if i != 1
           nextOlli!(basis_vector)
        end

        #ind1 = find_index(basis_vector)
        # this is literally the same as i

        for site in 1:L - 1
            modified_vector = copy(basis_vector)
            if basis_vector[site] > 0
                #HD[i, i] += basis_vector[site] * dis[site]
                modified_vector[site] -= 1
                modified_vector[site + 1] += 1
                index = find_indexOlli(modified_vector)
                HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1])
                HJ[index, i] = HJ[i, index]
                if basis_vector[site] > 1
                    HU[i, i] += -0.5 * basis_vector[site] * (basis_vector[site] - 1)
                end
            end
        end

        #HD[i, i] += basis_vector[L] * dis[L]
        if basis_vector[L] > 1
            HU[i, i] += -0.5 * basis_vector[L] * (basis_vector[L] - 1)
        end

        modified_vector = copy(basis_vector)
        if periodic == true && basis_vector[L] > 0
            modified_vector[L] -= 1
            modified_vector[1] += 1
            index = find_indexOlli(modified_vector)
            HJ[i, index] = sqrt(basis_vector[L] * modified_vector[1])
            HJ[index, i] = HJ[i, index]
        end
    end

    return HU, HJ
end


function interaction(L, N, dim; cap=N, dis = ones(L))
    HU = spzeros(dim, dim)
    basis_vector = first_state(L, N, cap)
    for i in 1:dim
        if i != 1
           next!(basis_vector, cap)
        end

        for site in 1:L
            if basis_vector[site] > 1
                HU[i, i] += -0.5 * basis_vector[site] * (basis_vector[site] - 1) * dis[L]
            end
        end

    end

    return HU
end


function hopping(L, N, dim; cap=N, periodic = false, dis = ones(L))
    HJ = spzeros(dim, dim)
    basis_vector = first_state(L, N, cap)
    for i in 1:dim
        if i != 1
           next!(basis_vector, cap)
        end

        for site in 1:L - 1
            modified_vector = copy(basis_vector)
            if basis_vector[site] > 0 && modified_vector[site + 1] + 1 <= cap
                modified_vector[site] -= 1
                modified_vector[site + 1] += 1
                index = find_index(modified_vector, cap)
                HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1]) * dis[L]
                HJ[index, i] = HJ[i, index]
            end
        end

        if periodic == true && basis_vector[L] > 0
            modified_vector = copy(basis_vector)
            modified_vector[L] -= 1
            modified_vector[1] += 1
            index = find_index(modified_vector, cap)
            HJ[i, index] = sqrt(basis_vector[L] * modified_vector[1])
            HJ[index, i] = HJ[i, index]
        end
    end

    return HJ
end


function disorder(L, N, dim; cap=N, return_min = false, dis = 2. .* rand(L) .- 1.)
    HD = spzeros(dim, dim)
    basis_vector = first_state(L, N, cap)
    for i in 1:dim
        if i != 1
           next!(basis_vector, cap)
        end

        for site in 1:L
            if basis_vector[site] > 0
                HD[i, i] += basis_vector[site] * dis[site]
            end
        end
    end

    if return_min == false
        return HD
    else
        return HD, argmin(dis)
    end
end


function reflect(L, N)
    dim = dimension(L, N)
    basis_vector = zeros(Int64, L)
    basis_vector[1] = N
    R = spzeros(dim, dim)
    for i in 1:dim
        if i != 1
           next!(basis_vector)
        end

        R[i, find_index(reverse(basis_vector))] = 1.
    end

    return R
end



function split_hamiltonian_2D(L1, L2, N; periodic = 0)
    sites = L1 * L2
    dim = dimension(sites, N)
    #HD = spzeros(dim, dim)
    HJ = spzeros(dim, dim)
    HU = spzeros(dim, dim)
    basis_vector = zeros(Int64, sites)
    basis_vector[1] = N

    if periodic == 0
        for i in 1:dim
            if i != 1
               next!(basis_vector)
            end

            for k in 1:L2
                for j in 1:L1
                    site = j + (k - 1) * L1
                    if basis_vector[site] > 0
                        if j < L1
                            modified_vector = copy(basis_vector)
                            modified_vector[site] -= 1
                            modified_vector[site + 1] += 1
                            index = find_index(modified_vector)
                            HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1])
                            HJ[index, i] = HJ[i, index]
                        end

                        if k < L2
                            modified_vector = copy(basis_vector)
                            modified_vector[site] -= 1
                            modified_vector[site + L1] += 1
                            index = find_index(modified_vector)
                            HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + L1])
                            HJ[index, i] = HJ[i, index]
                        end
                    end

                    if basis_vector[site] > 1
                        HU[i, i] += -0.5 * basis_vector[site] * (basis_vector[site] - 1)
                    end
                end
            end
        end

    elseif periodic == 1
        for i in 1:dim
            if i != 1
               next!(basis_vector)
            end

            for k in 1:L2
                for j in 1:L1
                    site = j + (k - 1) * L1
                    if basis_vector[site] > 0
                        if j < L1
                            modified_vector = copy(basis_vector)
                            modified_vector[site] -= 1
                            modified_vector[site + 1] += 1
                            index = find_index(modified_vector)
                            HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1])
                            HJ[index, i] = HJ[i, index]
                        else
                            modified_vector = copy(basis_vector)
                            modified_vector[site] -= 1
                            modified_vector[site + 1 - L1] += 1
                            index = find_index(modified_vector)
                            HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1 - L1])
                            HJ[index, i] = HJ[i, index]
                        end

                        if k < L2
                            modified_vector = copy(basis_vector)
                            modified_vector[site] -= 1
                            modified_vector[site + L1] += 1
                            index = find_index(modified_vector)
                            HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + L1])
                            HJ[index, i] = HJ[i, index]
                        end
                    end

                    if basis_vector[site] > 1
                        HU[i, i] += -0.5 * basis_vector[site] * (basis_vector[site] - 1)
                    end
                end
            end
        end
    elseif periodic == 2
        for i in 1:dim
            if i != 1
               next!(basis_vector)
            end

            for k in 1:L2
                for j in 1:L1
                    site = j + (k - 1) * L1
                    if basis_vector[site] > 0
                        if j < L1
                            modified_vector = copy(basis_vector)
                            modified_vector[site] -= 1
                            modified_vector[site + 1] += 1
                            index = find_index(modified_vector)
                            HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1])
                            HJ[index, i] = HJ[i, index]
                        else
                            modified_vector = copy(basis_vector)
                            modified_vector[site] -= 1
                            modified_vector[site + 1 - L1] += 1
                            index = find_index(modified_vector)
                            HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1 - L1])
                            HJ[index, i] = HJ[i, index]
                        end

                        if k < L2
                            modified_vector = copy(basis_vector)
                            modified_vector[site] -= 1
                            modified_vector[site + L1] += 1
                            index = find_index(modified_vector)
                            HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + L1])
                            HJ[index, i] = HJ[i, index]
                        else
                            modified_vector = copy(basis_vector)
                            modified_vector[site] -= 1
                            modified_vector[j] += 1
                            index = find_index(modified_vector)
                            HJ[i, index] = sqrt(basis_vector[site] * modified_vector[j])
                            HJ[index, i] = HJ[i, index]
                        end
                    end

                    if basis_vector[site] > 1
                        HU[i, i] += -0.5 * basis_vector[site] * (basis_vector[site] - 1)
                    end
                end
            end
        end
    end

    return HU, HJ
end


#~ function diag_hop(L, N)
#~     dim = dimension(L, N)
#~     #HD = spzeros(dim, dim)
#~     HJ = spzeros(dim, dim)
#~     HU = spzeros(dim, dim)
#~     basis_vector = zeros(Int64, L)
#~     basis_vector[1] = N
#~     Lp1 = 1.  / (1 + L)
#~     for i in 1:dim
#~         if i != 1
#~            next!(basis_vector)
#~         end

#~         for j in 1:L
#~             HJ[i, i] += 2. * cos(pi * j * Lp1) * basis_vector[j]
#~         end
#~     end

#~     return HJ
#~ end


#~ function rec_int_factor(L, k, l, m, n)
#~     out = 0
#~     if mod(k + l + m + n, L + 1) == 0
#~         out += (1 + (-1)^((k + l + m + n) / (L + 1)))
#~     end

#~     if mod(k + l - m - n, L + 1) == 0
#~         out += (1 + (-1)^((k + l - m - n) / (L + 1)))
#~     end

#~     if mod(k - l + m - n, L + 1) == 0
#~         out += (1 + (-1)^((k - l + m - n) / (L + 1)))
#~     end

#~     if mod(k - l - m + n, L + 1) == 0
#~         out += (1 + (-1)^((k - l - m + n) / (L + 1)))
#~     end

#~     if mod(k + l + m - n, L + 1) == 0
#~         out -= (1 + (-1)^((k + l + m - n) / (L + 1)))
#~     end

#~     if mod(k + l - m + n, L + 1) == 0
#~         out -= (1 + (-1)^((k + l - m + n) / (L + 1)))
#~     end

#~     if mod(k - l + m + n, L + 1) == 0
#~         out -= (1 + (-1)^((k - l + m + n) / (L + 1)))
#~     end

#~     if mod(k - l - m - n, L + 1) == 0
#~         out -= (1 + (-1)^((k - l - m - n) / (L + 1)))
#~     end

#~     return out / 2
#~ end


#~ function split_reciprocal_hamiltonian(L, N)
#~     dim = dimension(L, N)
#~     HJ = spzeros(dim, dim)
#~     HU = spzeros(dim, dim)
#~     basis_vector = zeros(Int64, L)
#~     basis_vector[1] = N
#~     for i in 1:dim
#~         if i != 1
#~            next!(basis_vector)
#~         end

#~         for k in 1:L
#~             HJ[i, i] += basis_vector[k] * cos(pi * k / (L + 1))
#~             if basis_vector[k] > 0
#~                 vector1 = copy(basis_vector)
#~                 n1 = vector1[k]
#~                 vector1[k] -= 1
#~                 for l in 1:L
#~                     if vector1[l] > 0
#~                         n2 = vector1[l]
#~                         vector1[l] -= 1
#~                         for m in 1:L
#~                             for n in 1:L
#~                                 vector2 = copy(vector1)
#~                                 if mod(k + l + m + n, 2) == 0
#~                                     vector2[m] += 1
#~                                     n3 = vector2[m]
#~                                     vector2[n] += 1
#~                                     n4 = vector2[n]

#~                                     j = find_index(vector2)
#~                                     fac = rec_int_factor(L, k, l, m, n)
#~                                     fac *= sqrt(n1 * n2 * n3 * n4)
#~                                     HU[i, j] = fac
#~                                 end
#~                             end
#~                         end
#~                     end
#~                 end
#~             end
#~         end
#~     end

#~     HJ *= 2.
#~     HU /= -2 * (L + 1)

#~     return HU, HJ
#~ end


#~ function hamiltonian(L, N, J, D)
#~     dim = dimension(L, N)
#~     H = spzeros(dim, dim)
#~     basis_vector = zeros(Int64, L)
#~     basis_vector[1] = N
#~     dis = 2. * rand(L) .- 1.
#~     for i in 1:dim
#~         if i != 1
#~            next!(basis_vector)
#~         end

#~         for site in 1:L - 1
#~             modified_vector = copy(basis_vector)
#~             if basis_vector[site] > 0
#~                 H[i, i] += D * basis_vector[site] * dis[site]
#~                 modified_vector[site] -= 1
#~                 modified_vector[site + 1] += 1
#~                 index = find_index(modified_vector)
#~                 H[i, index] = J * sqrt(basis_vector[site] * modified_vector[site + 1])
#~                 H[index, i] = H[i, index]
#~                 if basis_vector[site] > 1
#~                     H[i, i] += -0.5 * basis_vector[site] * (basis_vector[site] - 1)
#~                 end
#~             end
#~         end

#~         H[i, i] += basis_vector[L] * dis[L]
#~         if basis_vector[L] > 1
#~             H[i, i] += -0.5 * basis_vector[L] * (basis_vector[L] - 1)
#~         end
#~     end

#~     return H
#~ end




#~ function split_extended_hamiltonian(L, N)
#~     dim = dimension(L, N)
#~     #HD = spzeros(dim, dim)
#~     HJ = spzeros(dim, dim)
#~     HU = spzeros(dim, dim)
#~     basis_vector = zeros(Int64, L)
#~     basis_vector[1] = N
#~     for i in 1:dim
#~         if i != 1
#~            next!(basis_vector)
#~         end

#~         for site in 1:L - 1
#~             if basis_vector[site] > 0
#~                 modified_vector = copy(basis_vector)
#~                 #HD[i, i] += basis_vector[site] * dis[site]
#~                 modified_vector[site] -= 1
#~                 modified_vector[site + 1] += 1
#~                 index = find_index(modified_vector)
#~                 HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1])
#~                 HJ[index, i] = HJ[i, index]
#~                 if site <= L - 2
#~                     modified_vector = copy(basis_vector)
#~                     modified_vector[site] -= 1
#~                     modified_vector[site + 2] += 1
#~                     index = find_index(modified_vector)
#~                     HJ[i, index] = 0.1 * sqrt(basis_vector[site] * modified_vector[site + 2])
#~                     HJ[index, i] = HJ[i, index]
#~                 end

#~                 if basis_vector[site] > 1
#~                     HU[i, i] += -0.5 * basis_vector[site] * (basis_vector[site] - 1)
#~                     if basis_vector[site] > 2
#~                         HU[i, i] += basis_vector[site] * (basis_vector[site] - 1) * (basis_vector[site] - 2) / 60.0
#~                     end
#~                 end
#~             end
#~         end

#~         #HD[i, i] += basis_vector[L] * dis[L]
#~         if basis_vector[L] > 1
#~             HU[i, i] += -0.5 * basis_vector[L] * (basis_vector[L] - 1)
#~             if basis_vector[L] > 2
#~                 HU[i, i] += basis_vector[L] * (basis_vector[L] - 1) * (basis_vector[L] - 2) / 60.0
#~             end
#~         end
#~     end

#~     return HU, HJ
#~ end


#~ function split_extended_periodic_hamiltonian(L, N)
#~     dim = dimension(L, N)
#~     #HD = spzeros(dim, dim)
#~     HJ = spzeros(dim, dim)
#~     HU = spzeros(dim, dim)
#~     basis_vector = zeros(Int64, L)
#~     basis_vector[1] = N
#~     for i in 1:dim
#~         if i != 1
#~            next!(basis_vector)
#~         end

#~         for site in 1:L
#~             if basis_vector[site] > 1
#~                 HU[i, i] += -0.5 * basis_vector[site] * (basis_vector[site] - 1)
#~                 if basis_vector[site] > 2
#~                     HU[i, i] += basis_vector[site] * (basis_vector[site] - 1) * (basis_vector[site] - 2) / 60.0
#~                 end
#~             end

#~             if basis_vector[site] > 0
#~                 if site <= L - 2
#~                     modified_vector = copy(basis_vector)
#~                     #HD[i, i] += basis_vector[site] * dis[site]
#~                     modified_vector[site] -= 1
#~                     modified_vector[site + 1] += 1
#~                     index = find_index(modified_vector)
#~                     HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1])
#~                     HJ[index, i] = HJ[i, index]

#~                     modified_vector = copy(basis_vector)
#~                     modified_vector[site] -= 1
#~                     modified_vector[site + 2] += 1
#~                     index = find_index(modified_vector)
#~                     HJ[i, index] = 0.1 * sqrt(basis_vector[site] * modified_vector[site + 2])
#~                     HJ[index, i] = HJ[i, index]
#~                 elseif site == L - 1
#~                     modified_vector = copy(basis_vector)
#~                     #HD[i, i] += basis_vector[site] * dis[site]
#~                     modified_vector[site] -= 1
#~                     modified_vector[site + 1] += 1
#~                     index = find_index(modified_vector)
#~                     HJ[i, index] = sqrt(basis_vector[site] * modified_vector[site + 1])
#~                     HJ[index, i] = HJ[i, index]

#~                     modified_vector = copy(basis_vector)
#~                     modified_vector[site] -= 1
#~                     modified_vector[1] += 1
#~                     index = find_index(modified_vector)
#~                     HJ[i, index] = 0.1 * sqrt(basis_vector[site] * modified_vector[1])
#~                     HJ[index, i] = HJ[i, index]
#~                 elseif site == L
#~                     modified_vector = copy(basis_vector)
#~                     #HD[i, i] += basis_vector[site] * dis[site]
#~                     modified_vector[site] -= 1
#~                     modified_vector[1] += 1
#~                     index = find_index(modified_vector)
#~                     HJ[i, index] = sqrt(basis_vector[site] * modified_vector[1])
#~                     HJ[index, i] = HJ[i, index]

#~                     modified_vector = copy(basis_vector)
#~                     modified_vector[site] -= 1
#~                     modified_vector[2] += 1
#~                     index = find_index(modified_vector)
#~                     HJ[i, index] = 0.1 * sqrt(basis_vector[site] * modified_vector[2])
#~                     HJ[index, i] = HJ[i, index]
#~                 end
#~             end
#~         end
#~     end

#~     return HU, HJ
#~ end

#~ ----------------------------------------------- open system operators ---------------------------------------------------


function disorder_open(L, N; return_min = false, dis = 2. .* rand(L) .- 1.)
    dim = dimension_open(L, N)
    HD = spzeros(dim, dim)
    done = 1
    for n in 1:N
        ndim = dimension(L, n)
        basis_vector = zeros(Int64, L)
        basis_vector[1] = n
        for i in 1:ndim
            if i != 1
               next!(basis_vector)
            end

            for site in 1:L
                if basis_vector[site] > 0
                    HD[i + done, i + done] += basis_vector[site] * dis[site]
                end
            end
        end

        if return_min == false
            return HD
        else
            return HD, argmin(dis)
        end

        done += ndim
    end
end


function split_hamiltonian_open(L, N; periodic = false)
    dim = dimension_open(L, N)
    #HD = spzeros(dim, dim)
    HJ = spzeros(dim, dim)
    HU = spzeros(dim, dim)
    done = 1
    for n in 1:N
        ndim = dimension(L, n)
        basis_vector = zeros(Int64, L)
        basis_vector[1] = n
        for i in 1:ndim
            if i != 1
               next!(basis_vector)
            end

            for site in 1:L - 1
                modified_vector = copy(basis_vector)
                if basis_vector[site] > 0
                    modified_vector[site] -= 1
                    modified_vector[site + 1] += 1
                    index = find_index_open(modified_vector)
                    HJ[i + done, index] = sqrt(basis_vector[site] * modified_vector[site + 1])
                    HJ[index, i + done] = HJ[i + done, index]
                    if basis_vector[site] > 1
                        HU[i + done, i + done] += -0.5 * basis_vector[site] * (basis_vector[site] - 1)
                    end
                end
            end


            if basis_vector[L] > 1
                HU[i + done, i + done] += -0.5 * basis_vector[L] * (basis_vector[L] - 1)
            end

            if periodic == true && basis_vector[L] > 0
                modified_vector[L] -= 1
                modified_vector[1] += 1
                index = find_index_open(modified_vector)
                HJ[i + done, index] = sqrt(basis_vector[L] * modified_vector[1])
                HJ[index, i + done] = HJ[i + done, index]
            end
        end

        done += ndim
    end

    return HU, HJ
end


function annihilation(L, N)
    dim = dimension_open(L, N)
    a = spzeros(dim, dim)
    done = 1
    for n in 1:N
        ndim = dimension(L, n)
        basis_vector = zeros(Int64, L)
        basis_vector[1] = n
        for i in 1:ndim
            if i != 1
               next!(basis_vector)
            end

            for site in 1:L - 1
                modified_vector = copy(basis_vector)
                if basis_vector[site] > 0
                    modified_vector[site] -= 1
                    index = find_index_open(modified_vector)
                    a[index, i + done] = sqrt(basis_vector[site])
                end
            end
        end

        done += ndim
    end

    return a
end

function number_open(L, N, site)
    dim = dimension_open(L, N)
    n_l = spzeros(dim, dim)
    done = 1
    for n in 1:N
        ndim = dimension(L, n)
        basis_vector = zeros(Int64, L)
        basis_vector[1] = n
        for i in 1:ndim
            if i != 1
               next!(basis_vector)
            end

            n_l[i + done, i + done] += basis_vector[site]
        end

        done += ndim
    end

    return n_l
end


function numbers_open(L, N)
    n_l::Array{SparseMatrixCSC{Float64, Int64}, 1} = []
    for l in 1:L
        push!(n_l, number_open(L, N, l))
    end

    return n_l
end

"""
This module gives general methods for 2D momentum-block-diagonalized ED calculations.
Sectors of other quantum Numbers should be handled outside this module.
This module only sets sectors of total (crystal) momentum, also called blocks.
"""
module MomentumConservedExactDiagonalization

public MBS64, Scattering, ED_mbslist_onecomponent
export EDPara, ED_mbslist
export ED_momentum_block_division
export ED_sortedScatteringList_onebody
export ED_sortedScatteringList_twobody
export ED_Hmlt_Matrix, ED_Hmlt_Matrix_threaded, EDsolve

using LinearAlgebra
using KrylovKit
using ExtendableSparse, SparseArrays

# Include utilities
include("init_parameter.jl")
include("MBS.jl")
include("momentum_decomposition.jl")
include("scattering.jl")
include("sparse_matrix.jl")

















"""
Generate sorted lists of one-body scattering terms from the parameters.
Extracts one-body terms from EDpara.H1 for multi-component systems.
"""
function ED_sortedScatteringList_onebody(para::EDPara)
    sct_list1 = Vector{Scattering{1}}()
    Nk = para.Nk
    Nch = para.Nc_hopping
    Ncc = para.Nc_conserve

    # Extract one-body terms from H1[ch1, ch2, cc, k]
    for ch1 in 1:Nch, ch2 in 1:Nch, cc in 1:Ncc, k in 1:Nk
        V = para.H_onebody[ch1, ch2, cc, k]
        if !iszero(V)
            # Map component indices to global orbital indices
            i_ot = k + Nk * (ch1 - 1) + Nk * Nch * (cc - 1)  # output orbital
            i_in = k + Nk * (ch2 - 1) + Nk * Nch * (cc - 1)  # input orbital

            # Create scattering term with normal ordering
            i_in >= i_ot && push!(sct_list1, NormalScattering(V, i_ot, i_in))
        end
    end
    
    return sortMergeScatteringList(sct_list1)
end



"""
Generate grouped momentum pairs by their total momentum.
Uses indices to represent momenta instead of actual momentum values.
Output:
- Dict{Tuple{Int64,Int64}, Vector{Tuple{Int64,Int64}}}: 
Dictionary mapping total momentum to list of momentum index pairs (i,j)
"""
function group_momentum_pairs(para::EDPara)
    
    # Dictionary to store momentum groups
    momentum_groups = Dict{Tuple{Int64,Int64}, Vector{Tuple{Int64,Int64}}}()
    
    # Generate all possible pairs (including identical pairs)
    for i in 1:para.Nk, j in 1:i  # i >= j to avoid duplicates
        # Calculate total momentum using existing function
        K_total = MBS64_totalmomentum(para, i, j)
        pair_indices = (i, j)
        
        # Add to appropriate group
        if haskey(momentum_groups, K_total)
            push!(momentum_groups[K_total], pair_indices)
        else
            momentum_groups[K_total] = [pair_indices]
        end
    end
    
    return momentum_groups
end

"""
Generate all scattering terms between momentum pairs with the same total momentum.
According to EDplan.md step 2, for each pair of input/output momentum pairs,
we need to generate all 4 scattering terms accounting for state orderings.
"""
function scat_pair_group(pair_group::Vector{Tuple{Int64,Int64}}, para::EDPara;
    output::Bool = false)::Vector{Scattering{2}}
    
    scattering_list = Vector{Scattering{2}}()
    Nc = para.Nc
    Nk = para.Nk
    
    # Iterate over all input and output pairs
    for (ki1, ki2) in pair_group, (kf1, kf2) in pair_group
        output && println()
        output && println("ki1, ki2, kf1, kf2 = ($ki1, $ki2), ($kf1, $kf2)")
        # Generate all component index combinations
        for ci1 in 1:Nc, ci2 in 1:Nc, cf1 in 1:Nc, cf2 in 1:Nc
            
            # Map to global orbital indices
            # Global index = momentum_index + Nk * (component_index - 1)
            f1 = kf1 + Nk * (cf1 - 1)
            f2 = kf2 + Nk * (cf2 - 1)
            i1 = ki1 + Nk * (ci1 - 1)
            i2 = ki2 + Nk * (ci2 - 1)

            # no duplicate indices
            if i1 == i2 || f1 == f2
                continue
            end

            if ki1 == ki2 && i1 < i2
                continue
            end

            if kf1 == kf2 && f1 < f2
                continue
            end

            # inverse scattering only need to count onece, as the Hamiltonian is generated with upper half Hermitian()
            if minmax(i1, i2) >= minmax(f1, f2)

                # Calculate the direct and exchange amplitudes
                output && print(fldmod1(i1, Nk), fldmod1(i2, Nk), fldmod1(f1, Nk), fldmod1(f2, Nk),"        ")
                amp_direct = int_amp(i1, i2, f1, f2, para)
                output && print(fldmod1(i2, Nk), fldmod1(i1, Nk), fldmod1(f1, Nk), fldmod1(f2, Nk))
                amp_exchange = int_amp(i2, i1, f1, f2, para)
                amp = amp_direct - amp_exchange
                iszero(amp) || push!(scattering_list, NormalScattering(amp, f1, f2, i2, i1))
                output && println()
            end
        
        end
    end
    
    return scattering_list
end

"""
Generate sorted lists of two-body scattering terms from the parameters.
Uses V and F functions from EDPara to calculate scattering amplitudes.
This implements Step 3 from EDplan.md.
"""
function ED_sortedScatteringList_twobody(para::EDPara)

    sct_list2 = Vector{Scattering{2}}()
    
    momentum_groups = group_momentum_pairs(para)
    
    for (K_total, pairs) in momentum_groups
        append!(sct_list2, scat_pair_group(pairs, para))
    end
    
    return sortMergeScatteringList(sct_list2)
end







"""
Solve the Hamiltonian for the lowest n eigenvalues and eigenvectors

deal with possible errors.
"""
function matrix_solve(H::ExtendableSparseMatrix{ComplexF64}, n::Int64=6)
    vec0 = rand(ComplexF64, H.cscmatrix.m)
    vals, vecs, info = eigsolve(H, vec0, n, :SR, ishermitian=true)
    return vals, vecs
end

"""
warp function
"""
function EDsolve(sorted_mbs_block_list::Vector{MBS64{bits}}, 
    sorted_onebody_scat_list::Vector{Scattering{1}},
    sorted_twobody_scat_list::Vector{Scattering{2}}, n::Int64=6
) where {bits}
    # Construct H matrix
    H = ED_Hmlt_Matrix(sorted_mbs_block_list, 
        sorted_onebody_scat_list, sorted_twobody_scat_list)

    # solve H matrix
    vals, vecs = matrix_solve(H, n)

    # return results (energies, vectors)
    return vals, vecs
end




end

#=
"""
Calculate the reduced density matrix for a subsystem
"""
function reduced_density_matrix(psi::Vector{ComplexF64}, 
                               block_MBSList::Vector{MBS{bits}},
                               nA::Vector{Int64},
                               iA::Vector{Int64};
                               cutoff::Float64=1E-7) where {bits}
    
    bits_total = bits
    sorted_state_num_list = getfield.(block_MBSList, :n)
    
    # Filter by cutoff
    index_list = findall(x -> abs2(x) > cutoff, psi)
    psi_filtered = psi[index_list]
    num_list = sorted_state_num_list[index_list]
    
    Amask = sum(1 << i for i in iA; init=0)
    Bmask = (1 << bits_total) - 1 - Amask
    
    # Sort by B subsystem
    myless_fine(n1, n2) = n1 & Bmask < n2 & Bmask || n1 & Bmask == n2 & Bmask && n1 < n2
    myless_coarse(n1, n2) = n1 & Bmask < n2 & Bmask
    
    perm = sortperm(num_list; lt=myless_fine)
    psi_sorted = psi_filtered[perm]
    num_sorted = num_list[perm]
    
    # Find B chunks
    Bchunks_lastindices = Int64[0]
    let i=0
        while i < length(num_sorted)
            i = searchsortedlast(num_sorted, num_sorted[i+1]; lt=myless_coarse)
            push!(Bchunks_lastindices, i)
        end
    end
    
    NA = length(nA)
    RDM_threads = zeros(ComplexF64, NA, NA, Threads.nthreads())
    
    Threads.@threads for nchunk in 1:length(Bchunks_lastindices)-1
        id = Threads.threadid()
        chunkpiece = Bchunks_lastindices[nchunk]+1:Bchunks_lastindices[nchunk+1]
        numB = num_sorted[chunkpiece[1]] & Bmask
        
        for i in 1:NA
            numA = nA[i]
            num = numB + numA
            index = my_searchsortedfirst(num_sorted[chunkpiece], num)
            index == 0 && continue
            
            RDM_threads[i, i, id] += abs2(psi_sorted[chunkpiece[index]])
            
            for i′ in i+1:NA
                numA′ = nA[i′]
                num′ = numB + numA′
                index′ = my_searchsortedfirst(num_sorted[chunkpiece], num′)
                index′ == 0 && continue
                
                rhoii′ = conj(psi_sorted[chunkpiece[index′]]) * psi_sorted[chunkpiece[index]]
                RDM_threads[i, i′, id] += rhoii′
                RDM_threads[i′, i, id] += conj(rhoii′)
            end
        end
    end
    
    return sum(RDM_threads; dims=3)[:, :, 1]
end

"""
Calculate entanglement entropy from reduced density matrix
"""
function entanglement_entropy(RDM_A::Matrix{ComplexF64}; cutoff::Float64=1E-6)
    vals = eigvals(Hermitian(RDM_A))
    return sum(vals) do x
        if abs(x) < cutoff || x < 0
            return 0.0
        end
        -x * log2(x)
    end
end

"""
Calculate one-body reduced density matrix
"""
function one_body_reduced_density_matrix(psi::Vector{ComplexF64},
                                        block_MBSList::Vector{MBS{bits}};
                                        cutoff::Float64=1E-7) where {bits}
    
    bits_total = bits
    sorted_state_num_list = getfield.(block_MBSList, :n)
    
    # Filter by cutoff
    index_list = findall(x -> abs2(x) > cutoff, psi)
    psi_filtered = psi[index_list]
    num_list = sorted_state_num_list[index_list]
    
    # One-body RDM: ρ[i,j] = <c†_j c_i>
    rdm = zeros(ComplexF64, bits_total, bits_total)
    
    Threads.@threads for i in 0:bits_total-1
        for j in 0:bits_total-1
            if i == j
                # Diagonal: <n_i>
                for (idx, num) in enumerate(num_list)
                    if isodd(num >>> i)
                        rdm[i+1, j+1] += abs2(psi_filtered[idx])
                    end
                end
            else
                # Off-diagonal: <c†_j c_i>
                for (idx, num) in enumerate(num_list)
                    if isodd(num >>> i) && iseven(num >>> j)
                        new_num = num - (1 << i) + (1 << j)
                        new_idx = my_searchsortedfirst(num_list, new_num)
                        new_idx == 0 && continue
                        
                        sign_flip = sum(k -> (num >>> k) % 2, min(i,j)+1:max(i,j)-1; init=0)
                        sign = isodd(sign_flip) ? -1 : 1
                        
                        rho_ij = sign * conj(psi_filtered[new_idx]) * psi_filtered[idx]
                        rdm[i+1, j+1] += rho_ij
                        rdm[j+1, i+1] += conj(rho_ij)
                    end
                end
            end
        end
    end
    
    return rdm
end

# Export all public functions
public EDPara, ED_mbslist, ED_totalmomentum, ED_momentum_block_division
public EDHamiltonian, EDsolve, add_one_body!, add_interaction!
public reduced_density_matrix, entanglement_entropy, one_body_reduced_density_matrix
public normalize_interaction_list, normal_order!


=#
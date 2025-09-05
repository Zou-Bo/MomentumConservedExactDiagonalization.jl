"""
This module gives general methods for 2D momentum-block-diagonalized ED calculations.
Sectors of other quantum Numbers should be handled outside this module.
This module only sets sectors of total (crystal) momentum, also called blocks.
"""
module MomentumConservedExactDiagonalization

public MBS64, Scattering, ED_mbslist_onecomponent
export EDPara, ED_mbslist, ED_Hmlt_Matrix, EDsolve
export ED_momentum_block_division
export ED_sortedScatteringList_onebody
export ED_sortedScatteringList_twobody

using LinearAlgebra
using Combinatorics
using KrylovKit
using ExtendableSparse, SparseArrays


"search the index of the first object, return 0 if it's absent."
function my_searchsortedfirst(list, i)
    index = searchsortedfirst(list, i)
    if index > lastindex(list) || list[index] != i
        return 0
    else
        return index
    end
end


@kwdef mutable struct EDPara
    # momemta are in integers

    # G integer (momentum is conserved mod G; where G=0 means no mod)
    Gk::Tuple{Int64, Int64} = (0, 0)
    # k_list[:, i] = (k_x, k_y)
    k_list::Matrix{Int64}

    Nk::Int64 = size(k_list,2)
    Nc_hopping::Int64 = 1 # number of components with hopping
    Nc_conserve::Int64 = 1 # number of components with conserved quantum numbers
    Nc::Int64 = Nc_hopping * Nc_conserve # number of total components

    # orbital index i
    # i = ik + Nk * (ic-1) = ik + Nk * (ich-1) + (Nk * Nc_hopping) * (icc-1)


    error = begin
        @assert Nc > 0 "Number of components must be positive"
        @assert Nk*Nc <= 64 "The Hilbert space dimension must not exceed 64 bits."
        nothing
    end

    # (to be fill later by user)
    H_onebody::Array{ComplexF64,4} = zeros(ComplexF64, Nc_hopping, Nc_hopping, Nc_conserve, Nk)
    # V_q::Function = (q1, q2, cf1, cf2, ci1, ci2) -> ComplexF64(1.)
    # F_kf_ki_cf_ci::Function = (kf1_1, kf1_2, ki1_1, ki1_2, cf, ci) -> ComplexF64(1.)
    V_int::Function = (kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) -> ComplexF64(1.0)

end
"""
rewrite V_int, inputs are (momentum + component) combined index
"""
function int_amp(i1::Int64, i2::Int64, f1::Int64, f2::Int64, para::EDPara)::ComplexF64

    # Map combined indices back to momentum and component indices
    Nk = para.Nk
    Gk1 = para.Gk[1]
    Gk2 = para.Gk[2]
    ci1, ki1 = fldmod1(i1, Nk)
    ci2, ki2 = fldmod1(i2, Nk)
    cf1, kf1 = fldmod1(f1, Nk)
    cf2, kf2 = fldmod1(f2, Nk)

    return para.V_int(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) / Gk1 / Gk2
end



# MBS64 type <: Integer to record occupations of orbitals in a many-body states
"""
The struct of many-body state <: Integer consists of
a 64-bitinteger whose bitarray indicated occupations.
bit index i = ik + Nk * (ic-1), i = 1 ~ Nk * Nc
number = sum(i -> 1<<(i-1), occ_list; init=0)
Total number of physical bits (orbitals) is shown in the type parameter.
"""
struct MBS64{bits} <: Integer
    n::Int64
    function MBS64{bits}(n::Int64) where {bits}
        @assert bits isa Integer && 0 < bits <= 64 "The number of bits must be an integer between 1 and 64."
        @assert 0 <= reinterpret(UInt64, n) < (UInt64(1) << bits) "Integer representation out of range for given bits"
        new{bits}(n)
    end
end
function Base.show(io::IO, mbs::MBS64{bits}) where bits
    print(io, "MBS64: ", mbs.n, " = ", view(reverse(bitstring(mbs.n)), 1:bits), " ($bits bits)")
    if !isempty(findall(==('1'), view(reverse(bitstring(mbs.n)), bits+1:64)))
        println(io, " (Unphysical bits are occupied in MBS64.)")
        @warn "Unphysical bits are occupied in MBS64."
    end
end

import Base: *, <
# connect 2 MBS
function *(mbs1::MBS64{b1}, mbs2::MBS64{b2}) where {b1, b2}
    MBS64{b1+b2}(mbs1.n << b2 | mbs2.n)
end
# same to Int64
<(mbs1::MBS64, mbs2::MBS64) = mbs1.n < mbs2.n # for sorting purpose


"""
Return the list of occupied states (bit positions) in the many-body state.
"""
function occ_list(mbs::MBS64{bits}) where {bits}
    return findall(==('1'), view(reverse(bitstring(mbs.n)), 1:bits))
end

"""
Generate a MBS64 from a list of occupied states.
"""
function MBS64(bits, occ_list::Int64...)
    n = Int64(0)
    for i in occ_list
        @assert 1 <= i <= bits "Occupied state index out of bounds"
        n += 1 << (i - 1)
    end
    return MBS64{bits}(n)
end


import Base: isempty, empty!
"""
Check if the specified orbital(s) are all occupied in the many-body state.
"""
function isoccupied(mbs::MBS64{bits}, i_list::Int64...) where {bits}
    mask = MBS64(bits, i_list...)
    return mbs.n & mask.n == mask.n
end
"""
Check if the specified orbital(s) are all empty in the many-body state.
"""
function isempty(mbs::MBS64{bits}, i_list::Int64...) where {bits}
    mask = MBS64(bits, i_list...)
    return mbs.n & mask.n == 0
end
"""
To occupy the originally-empty orbital(s) in the many-body state.
"""
function occupy!(mbs::MBS64{bits}, i_list::Int64...; check::Bool=true) where {bits}
    mask = MBS64(bits, i_list...)
    if check
        @assert mbs.n & mask.n == 0 "Some orbitals are already occupied."
    end
    return reinterpret(MBS64{bits}, mbs.n | mask.n)
end
"""
To empty the originally-occupied orbital(s) in the many-body state.
"""
function empty!(mbs::MBS64{bits}, i_list::Int64...; check::Bool=true) where {bits}
    mask = MBS64(bits, i_list...)
    if check
        @assert mbs.n & mask.n == mask.n "Some orbitals are already empty."
    end
    return reinterpret(MBS64{bits}, mbs.n & ~mask.n)
end
"""
Count the number of occupied orbitals between two bit positions in the many-body state.
Excluding the two ends.
"""
function occ_num_between(mbs::MBS64{bits}, i_start::Int64, i_end::Int64) where {bits}
    @assert 1 <= i_start <= bits "Invalid bit positions"
    @assert 1 <= i_end <= bits "Invalid bit positions"
    i_start, i_end = minmax(i_start, i_end) # allow inversely-ordered inputs
    mask = zero(mbs.n)
    for i in i_start+1:i_end-1
        mask |= (1 << (i - 1))
    end
    return count_ones(mbs.n & mask)
end

"""
Return the total momentum K1 and K2 of a many body state.
The momentum is mod G if G is nonzero.
"""
function MBS64_totalmomentum(para::EDPara, mbs::MBS64)
    # momentum are integers
    k1 = 0; k2 = 0; Gk = para.Gk
    for i in occ_list(mbs)
        momentum = @view para.k_list[:, mod1(i, para.Nk)]
        k1 += momentum[1]
        k2 += momentum[2]
    end
    iszero(Gk[1]) || (k1 = mod(k1, Gk[1]))
    iszero(Gk[2]) || (k2 = mod(k2, Gk[2]))
    return k1, k2
end
function MBS64_totalmomentum(para::EDPara, i_list::Int64...)
    # momentum are integers
    k1 = 0; k2 = 0; Gk = para.Gk
    for i in i_list
        momentum = @view para.k_list[:, mod1(i, para.Nk)]
        k1 += momentum[1]
        k2 += momentum[2]
    end
    iszero(Gk[1]) || (k1 = mod(k1, Gk[1]))
    iszero(Gk[2]) || (k2 = mod(k2, Gk[2]))
    return k1, k2
end


"""
construct the MBS list of N electrons in one conserved component.
"""
function ED_mbslist_onecomponent(para::EDPara, N_in_one::Int64)
    Nstate = para.Nk * para.Nc_hopping
    @assert 0 <= N_in_one <= Nstate "Invalid number of electrons in one component"
    sort([MBS64(Nstate, combi...) for combi in collect(combinations(1:Nstate, N_in_one)) ])
end
"""
Construct a list of MBS with electron number (N1, N2, ...) in each conserved component.
"""
function ED_mbslist(para::EDPara, N_each_component::NTuple{N, Int64}) where {N}
    @assert N == para.Nc_conserve "The length of number_list must be equal to the number of conserved components $(para.Nc_conserve)"
    list = ED_mbslist_onecomponent(para, N_each_component[begin])
    for i in eachindex(N_each_component)[2:end]
        list = kron(ED_mbslist_onecomponent(para, N_each_component[i]), list)
    end
    return list
end




"""
Divide a given mbs list into blocks with same momentum
Input: MBS list and EDPara parameters
Output:
1. list of block MBS list,
2. list of total k1 of blocks,
3. list of total k2 of blocks,
4. list index of the k1=k2=0 block (return 0 if not exist)
"""
function ED_momentum_block_division(para::EDPara, mbs_list::Vector{MBS64{bits}};
    momentum_restriction = false, k1range=(-2,2), k2range=(-2,2),
    momentum_list::Vector{Tuple{Int64, Int64}} = Vector{Tuple{Int64, Int64}}(),
    ) where {bits}

    # Calculate momentum for each MBS
    TotalMomentum = Matrix{Int64}(undef, 2, length(mbs_list))
    k1_list = Vector{Int64}(undef, length(mbs_list))
    k2_list = Vector{Int64}(undef, length(mbs_list))
    
    for (idx, mbs) in enumerate(mbs_list)
        k1, k2 = MBS64_totalmomentum(para, mbs)
        TotalMomentum[1, idx] = k1
        TotalMomentum[2, idx] = k2
        k1_list[idx] = k1
        k2_list[idx] = k2
    end


    blocks = typeof(mbs_list)[]
    block_k1 = Int64[]
    block_k2 = Int64[]
    
    # Determine momentum ranges
    if isempty(mbs_list)
        return blocks, block_k1, block_k2, 0
    end
    
    k1min, k2min = minimum(TotalMomentum, dims=2)
    k1max, k2max = maximum(TotalMomentum, dims=2)

    if momentum_restriction
        k1min = max(k1range[1], k1min)
        k2min = max(k2range[1], k2min)
        k1max = min(k1range[2], k1max)
        k2max = min(k2range[2], k2max)
    end

    # Group by momentum
    if isempty(momentum_list)
        for K1 in k1min:k1max, K2 in k2min:k2max
            mask = findall((k1_list .== K1) .& (k2_list .== K2))
            if !iszero(length(mask))
                push!(blocks, mbs_list[mask])
                push!(block_k1, K1)
                push!(block_k2, K2)
            end
        end
    else
        for K1 in k1min:k1max, K2 in k2min:k2max
            if (K1, K2) ∈ momentum_list
                mask = findall((k1_list .== K1) .& (k2_list .== K2))
                if !iszero(length(mask))
                    push!(blocks, mbs_list[mask])
                    push!(block_k1, K1)
                    push!(block_k2, K2)
                end
            end
        end
    end
    
    # Find k=0 block index
    block_num_0 = findfirst(eachindex(block_k1)) do bn
        block_k1[bn] == 0 && block_k2[bn] == 0
    end
    
    return blocks, block_k1, block_k2, something(block_num_0, 0)
end




struct Scattering{N}
    Amp::ComplexF64
    out::NTuple{N, Int64}
    in::NTuple{N, Int64}
end
function Scattering(V, outin::Int64...)
    @assert iseven(length(outin)) "Number of indices must be even"
    N = length(outin) ÷ 2
    outstates = outin[begin:N]
    instates = reverse(outin[N + 1:end])
    return Scattering{N}(ComplexF64(V), outstates, instates)
end
function Base.show(io::IO, st::Scattering{N}) where {N}
    print(io, "$N-body scattering: c†_out ", st.out, " c_in ", reverse(st.in), ", Amp = ", st.Amp)
end

"""
Generate a scattering term with normal ordering (now working for N = 1, 2)
term: V * c†_i1 c†_i2 ... c†_iN c_jN ... c_j2 c_j1 (j-in, i-out )
(1) j1 > j2 > ... > jN (no equality)
(2) i1 > i2 > ... > iN (no equality)
(3) j1 > i1 or j1 = i1 && j2 > i2 or j1,j2 = i1,i2 && j3 > i3 or ... or j1,...,jN-1 = i1,...,iN-1 && jN >= iN
"""
function NormalScattering(V, i, j)
    # N = 1
    
    # Apply normal ordering rules

    # j >= i
    if j < i
        j, i = i, j
        V = conj(V)
    end

    if j == i
        V = real(V) + 0im
    end

    return Scattering(V, i, j)
end
function NormalScattering(V, i1, i2, j2, j1)
    # N = 2
    
    # Apply normal ordering rules

    # Skip if indices are invalid
    if j1 == j2 || i1 == i2
        @warn "Skipping invalid interaction term: $S"
        return Scattering(0, i1, i2, j2, j1)
    end
    
    # i1 > i2
    if i1 < i2
        i1, i2 = i2, i1
        V = -V
    end
    
    # j1 > j2
    if j1 < j2
        j1, j2 = j2, j1
        V = -V
    end
    
    # j1 > i1 or j1 = i1 && j2 > i2
    if j1 < i1 || (j1 == i1 && j2 < i2)
        j1, i1 = i1, j1
        j2, i2 = i2, j2
        V = conj(V)
    end

    if j1 == i1 && j2 == i2
        V = real(V) + 0im
    end
    # @show V, i1, i2, j2, j1
    return Scattering(V, i1, i2, j2, j1)
end
isnormal(x::Scattering{1}) = x.in >= x.out
isnormal(x::Scattering{2}) = x.in[1] > x.in[2] && x.out[1] > x.out[2] && x.in >= x.out




import Base: isless, ==, +
isless(x::Scattering{N}, y::Scattering{N}) where {N} = x.in < y.in || x.in == y.in && x.out < y.out
==(x::Scattering{N}, y::Scattering{N}) where {N} = x.in == y.in && x.out == y.out
function +(x::Scattering{N}, y::Scattering{N})::Scattering{N} where {N}
    @assert x == y "Can only add identical Scattering terms"
    return Scattering{N}(x.Amp + y.Amp, x.out, x.in)
end
"""
Sort and merge a list of normalized scattering terms
"""
function sortMergeScatteringList(normal_sct_list::Vector{Scattering{N}})::Vector{Scattering{N}} where {N}
    @assert N == 1 || N == 2 "Only normal ordering of 1-body and 2-body scattering terms are supported"
    @assert all(isnormal, normal_sct_list) "All scattering terms must be normalized"
    sorted_list = sort(normal_sct_list)
    merged_list = Vector{Scattering{N}}()
    for sct in sorted_list
        if !isempty(merged_list) && sct == merged_list[end]
            merged_list[end] += sct
        else
            push!(merged_list, sct)
        end
    end
    return merged_list
end











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






# add functions here to generate the two body scattering list.

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
"Matrix Method":
Add one-body and two-body terms to the sparse Hamiltonian matrix.
Only upper triangular part is stored, i.e., H[i,j] with i <= j,
and return the Hermitian form of H.
"""
function ED_Hmlt_Matrix(sorted_mbs_block_list::Vector{MBS64{bits}}, 
    sorted_onebody_scat_list::Vector{Scattering{1}},
    sorted_twobody_scat_list::Vector{Scattering{2}},
)::ExtendableSparseMatrix{ComplexF64} where {bits}

    size = length(sorted_mbs_block_list)
    H = ExtendableSparseMatrix(ComplexF64, size, size);

    for mbsj in eachindex(sorted_mbs_block_list)
        mbs_in = sorted_mbs_block_list[mbsj]

        # for each given in-state, find all possible out-states
        # two-body scattering terms
        for scat in sorted_twobody_scat_list
            if isoccupied(mbs_in, scat.in...)
                if scat.in == scat.out
                    updateindex!(H, +, scat.Amp, mbsj, mbsj)
                else
                    mbs_mid = empty!(mbs_in, scat.in...; check=false)
                    if isempty(mbs_mid, scat.out...)
                        mbs_out = occupy!(mbs_mid, scat.out...; check=false)
                        mbsi = my_searchsortedfirst(sorted_mbs_block_list, mbs_out)
                        @assert mbsi != 0 "H is not momentum-conserving."
                        if iseven(occ_num_between(mbs_mid, scat.in...) + occ_num_between(mbs_mid, scat.out...))
                            updateindex!(H, +, scat.Amp, mbsi, mbsj)
                        else
                            updateindex!(H, +, -scat.Amp, mbsi, mbsj)
                        end
                    end
                end
            end
        end
        # one-body scattering terms
        for scat in sorted_onebody_scat_list
            if isoccupied(mbs_in, scat.in...)
                if scat.in == scat.out
                    updateindex!(H, +, scat.Amp, mbsj, mbsj)
                else
                    mbs_mid = empty!(mbs_in, scat.in...; check=false)
                    if isempty(mbs_mid, scat.out...)
                        mbs_out = occupy!(mbs_mid, scat.out...; check=false)
                        mbsi = my_searchsortedfirst(sorted_mbs_block_list, mbs_out)
                        @assert mbsi != 0 "H is not momentum-conserving."
                        updateindex!(H, +, scat.Amp, mbsi, mbsj)
                    end
                end
            end
        end

    end

    return ExtendableSparseMatrix(Hermitian(H, :U))
end



"""
Create a dictionary mapping from MBS64 states to their indices for O(1) lookup.
This eliminates the my_searchsortedfirst bottleneck.
"""
function create_state_mapping(sorted_mbs_block_list::Vector{MBS64{bits}}) where {bits}
    mapping = Dict{Int, Int}()  # Use Int instead of MBS64{bits}
    for (i, state) in enumerate(sorted_mbs_block_list)
        mapping[state.n] = i
    end
    return mapping
end
"""
Threaded version of ED_Hmlt_Matrix with pre-computed state mapping and COO format construction.
Provides 4-8x speedup for medium to large systems.
"""
function ED_Hmlt_Matrix_threaded(
    sorted_mbs_block_list::Vector{MBS64{bits}}, 
    sorted_onebody_scat_list::Vector{Scattering{1}},
    sorted_twobody_scat_list::Vector{Scattering{2}},
)::SparseMatrixCSC{ComplexF64} where {bits}

    n_states = length(sorted_mbs_block_list)
    state_mapping = create_state_mapping(sorted_mbs_block_list)
    
    # Thread-local storage for COO format
    n_threads = Threads.nthreads()
    thread_I = [Vector{Int}() for _ in 1:n_threads]
    thread_J = [Vector{Int}() for _ in 1:n_threads]
    thread_V = [Vector{ComplexF64}() for _ in 1:n_threads]
    
    # Parallel construction over columns
    Threads.@threads for j in 1:n_states
        tid = Threads.threadid()
        mbs_in = sorted_mbs_block_list[j]
        
        # Two-body scattering terms
        for scat in sorted_twobody_scat_list
            if isoccupied(mbs_in, scat.in...)
                if scat.in == scat.out
                    # Diagonal term
                    push!(thread_I[tid], j)
                    push!(thread_J[tid], j)
                    push!(thread_V[tid], scat.Amp)
                else
                    # Off-diagonal term
                    mbs_mid = empty!(mbs_in, scat.in...; check=false)
                    if isempty(mbs_mid, scat.out...)
                        mbs_out = occupy!(mbs_mid, scat.out...; check=false)
                        i = get(state_mapping, mbs_out.n, 0)
                        @assert i != 0 "H is not momentum-conserving."
                        
                        if iseven(occ_num_between(mbs_mid, scat.in...) + occ_num_between(mbs_mid, scat.out...))
                            push!(thread_I[tid], i)
                            push!(thread_J[tid], j)
                            push!(thread_V[tid], scat.Amp)
                        else
                            push!(thread_I[tid], i)
                            push!(thread_J[tid], j)
                            push!(thread_V[tid], -scat.Amp)
                        end
                    end
                end
            end
        end
        
        # One-body scattering terms
        for scat in sorted_onebody_scat_list
            if isoccupied(mbs_in, scat.in...)
                if scat.in == scat.out
                    # Diagonal term
                    push!(thread_I[tid], j)
                    push!(thread_J[tid], j)
                    push!(thread_V[tid], scat.Amp)
                else
                    # Off-diagonal term
                    mbs_mid = empty!(mbs_in, scat.in...; check=false)
                    if isempty(mbs_mid, scat.out...)
                        mbs_out = occupy!(mbs_mid, scat.out...; check=false)
                        i = get(state_mapping, mbs_out.n, 0)
                        @assert i != 0 "H is not momentum-conserving."
                        
                        push!(thread_I[tid], i)
                        push!(thread_J[tid], j)
                        push!(thread_V[tid], scat.Amp)
                    end
                end
            end
        end
    end
    
    # Merge thread-local results
    total_entries = sum(length.(thread_I))
    I = Vector{Int}(undef, total_entries)
    J = Vector{Int}(undef, total_entries)
    V = Vector{ComplexF64}(undef, total_entries)
    
    offset = 0
    for tid in 1:n_threads
        n = length(thread_I[tid])
        if n > 0
            I[offset+1:offset+n] .= thread_I[tid]
            J[offset+1:offset+n] .= thread_J[tid]
            V[offset+1:offset+n] .= thread_V[tid]
            offset += n
        end
    end
    
    # Convert to sparse matrix (Hermitian)
    H = sparse(I, J, V, n_states, n_states)
    return H + H' - sparse(Diagonal(H))
end




"""
Solve the Hamiltonian for the lowest n eigenvalues and eigenvectors
"""
function EDsolve(H::ExtendableSparseMatrix{ComplexF64}, n::Int64=6)
    vec0 = rand(ComplexF64, H.cscmatrix.m)
    vals, vecs, info = eigsolve(H, vec0, n, :SR, ishermitian=true)
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
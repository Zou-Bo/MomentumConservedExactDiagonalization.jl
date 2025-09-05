"""
    Hamiltonian.jl - Scattering-based Hamiltonian construction
    
    This module provides the scattering list generation system for momentum-conserved
    exact diagonalization calculations, including dedicated structs for 1-body and
    2-body scattering processes.
"""

export Scattering1Body, Scattering2Body
export generate_onebody_scattering, generate_twobody_scattering
export sort_scattering_list, filter_zero_amplitude

using LinearAlgebra

"""
    Scattering1Body

Represents a one-body scattering process with initial and final orbital indices
and scattering amplitude.

# Fields
- `i::Int64`: Initial orbital index
- `f::Int64`: Final orbital index  
- `amp::ComplexF64`: Scattering amplitude

# Orbital Indexing
Uses the formula: `i = ik + Nk * (ic-1) + (Nk * Nc_hopping) * (icc-1)`
where `ik` is momentum index, `ic` is component index, `icc` is conserved component index
"""
struct Scattering1Body
    i::Int64  # Initial orbital index
    f::Int64  # Final orbital index
    amp::ComplexF64  # Scattering amplitude
end

"""
    Scattering2Body

Represents a two-body scattering process with initial and final orbital indices
and scattering amplitude. Momentum conservation: `kf1 + kf2 = ki1 + ki2` (mod G)

# Fields
- `i1::Int64`: Initial orbital 1
- `i2::Int64`: Initial orbital 2
- `f1::Int64`: Final orbital 1
- `f2::Int64`: Final orbital 2
- `amp::ComplexF64`: Scattering amplitude

# Momentum Conservation
Ensures `kf1 + kf2 = ki1 + ki2` (mod G) where G is from EDPara.Gk
"""
struct Scattering2Body
    i1::Int64  # Initial orbital 1
    i2::Int64  # Initial orbital 2
    f1::Int64  # Final orbital 1
    f2::Int64  # Final orbital 2
    amp::ComplexF64  # Scattering amplitude
end

# Import EDPara for type stability
using ..MomentumConservedExactDiagonalization: EDPara

"""
    orbital_index(k::Int64, ch::Int64, cc::Int64, para::EDPara)::Int64

Calculate global orbital index from momentum and component indices.
Uses formula: `i = k + Nk * (ch-1) + (Nk * Nc_hopping) * (cc-1)`
"""
function orbital_index(k::Int64, ch::Int64, cc::Int64, para::EDPara)::Int64
    return k + para.Nk * (ch - 1) + (para.Nk * para.Nc_hopping) * (cc - 1)
end

"""
    decompose_orbital_index(i::Int64, para::EDPara)::Tuple{Int64, Int64, Int64}

Decompose global orbital index into momentum and component indices.
Returns (k, ch, cc) where:
- k: momentum index (1:Nk)
- ch: hopping component index (1:Nc_hopping)  
- cc: conserved component index (1:Nc_conserve)
"""
function decompose_orbital_index(i::Int64, para::EDPara)::Tuple{Int64, Int64, Int64}
    Nk = para.Nk
    Nch = para.Nc_hopping
    
    # First get conserved component
    cc, remainder = fldmod1(i, Nk * Nch)
    
    # Then get hopping component and momentum
    ch, k = fldmod1(remainder, Nk)
    
    return (k, ch, cc)
end

"""
    is_upper_triangle(i::Int64, f::Int64)::Bool

Check if orbital indices satisfy upper triangle condition (i <= f).
This ensures we only generate unique matrix elements for Hermitian matrices.
"""
function is_upper_triangle(i::Int64, f::Int64)::Bool
    return i <= f
end

"""
    is_upper_triangle_2body(i1::Int64, i2::Int64, f1::Int64, f2::Int64)::Bool

Check if two-body orbital indices satisfy upper triangle condition.
For 2-body terms, we use lexicographic ordering: (i1,i2) <= (f1,f2)
"""
function is_upper_triangle_2body(i1::Int64, i2::Int64, f1::Int64, f2::Int64)::Bool
    # Ensure proper ordering within each pair
    i_min, i_max = minmax(i1, i2)
    f_min, f_max = minmax(f1, f2)
    
    # Lexicographic comparison: (i_min, i_max) <= (f_min, f_max)
    return (i_min < f_min) || (i_min == f_min && i_max <= f_max)
end

"""
    apply_momentum_conservation(ki1::Int64, ki2::Int64, kf1::Int64, kf2::Int64, para::EDPara)::Bool

Check if momentum conservation is satisfied: `kf1 + kf2 = ki1 + ki2` (mod G)
"""
function apply_momentum_conservation(ki1::Int64, ki2::Int64, kf1::Int64, kf2::Int64, para::EDPara)::Bool
    Gk1, Gk2 = para.Gk
    
    # Calculate total initial and final momentum
    k1_total_init = para.k_list[1, ki1] + para.k_list[1, ki2]
    k2_total_init = para.k_list[2, ki1] + para.k_list[2, ki2]
    
    k1_total_final = para.k_list[1, kf1] + para.k_list[1, kf2]
    k2_total_final = para.k_list[2, kf1] + para.k_list[2, kf2]
    
    # Apply momentum conservation (mod G)
    if Gk1 != 0
        k1_total_init = mod(k1_total_init, Gk1)
        k1_total_final = mod(k1_total_final, Gk1)
    end
    
    if Gk2 != 0
        k2_total_init = mod(k2_total_init, Gk2)
        k2_total_final = mod(k2_total_final, Gk2)
    end
    
    return k1_total_init == k1_total_final && k2_total_init == k2_total_final
end

"""
    generate_onebody_scattering(para::EDPara)::Vector{Scattering1Body}

Generate one-body scattering list from H_onebody array in EDPara.
Extracts non-zero terms and applies upper triangle optimization.

# Process
1. Iterate through H_onebody[cf, ci, cc, k] array
2. Map component indices to global orbital indices
3. Apply upper triangle condition (i <= f)
4. Skip zero amplitude terms
5. Return sorted list by (i, f) for efficient processing
"""
function generate_onebody_scattering(para::EDPara)::Vector{Scattering1Body}
    scattering_list = Vector{Scattering1Body}()
    
    Nk = para.Nk
    Nch = para.Nc_hopping
    Ncc = para.Nc_conserve
    
    # Extract one-body terms from H_onebody[cf, ci, cc, k]
    for ch_f in 1:Nch, ch_i in 1:Nch, cc in 1:Ncc, k in 1:Nk
        amp = para.H_onebody[ch_f, ch_i, cc, k]
        
        # Skip zero amplitude terms
        if iszero(amp)
            continue
        end
        
        # Map component indices to global orbital indices
        i_orb = orbital_index(k, ch_i, cc, para)  # initial orbital
        f_orb = orbital_index(k, ch_f, cc, para)  # final orbital
        
        # Apply upper triangle optimization and component conservation
        if is_upper_triangle(i_orb, f_orb)
            push!(scattering_list, Scattering1Body(i_orb, f_orb, amp))
        end
    end
    
    # Sort by (i, f) for consistent ordering
    sort!(scattering_list, by = s -> (s.i, s.f))
    
    return scattering_list
end

"""
    generate_twobody_scattering(para::EDPara)::Vector{Scattering2Body}

Generate two-body scattering list using V_int function from EDPara.
Evaluates function for all momentum-conserving combinations with upper triangle optimization.

# Process
1. Generate all momentum pairs with same total momentum
2. For each pair combination, evaluate V_int function
3. Apply momentum conservation: kf1 + kf2 = ki1 + ki2 (mod G)
4. Apply component conservation for Nc_conserve quantum numbers
5. Use upper triangle optimization to avoid duplicates
6. Return sorted list for efficient processing
"""
function generate_twobody_scattering(para::EDPara)::Vector{Scattering2Body}
    scattering_list = Vector{Scattering2Body}()
    
    Nk = para.Nk
    Nc = para.Nc  # Total components = Nc_hopping * Nc_conserve
    
    # Group momentum pairs by total momentum for efficient processing
    momentum_groups = Dict{Tuple{Int64,Int64}, Vector{Tuple{Int64,Int64}}}()
    
    for ki1 in 1:Nk, ki2 in 1:Nk
        # Calculate total momentum
        k1_total = para.k_list[1, ki1] + para.k_list[1, ki2]
        k2_total = para.k_list[2, ki1] + para.k_list[2, ki2]
        
        # Apply momentum conservation (mod G)
        Gk1, Gk2 = para.Gk
        if Gk1 != 0
            k1_total = mod(k1_total, Gk1)
        end
        if Gk2 != 0
            k2_total = mod(k2_total, Gk2)
        end
        
        K_total = (k1_total, k2_total)
        
        if !haskey(momentum_groups, K_total)
            momentum_groups[K_total] = Vector{Tuple{Int64,Int64}}()
        end
        push!(momentum_groups[K_total], (ki1, ki2))
    end
    
    # Generate scattering terms for each momentum group
    for (K_total, momentum_pairs) in momentum_groups
        # For each input and output pair with same total momentum
        for (ki1, ki2) in momentum_pairs, (kf1, kf2) in momentum_pairs
            
            # Generate all component combinations
            for ci1 in 1:Nc, ci2 in 1:Nc, cf1 in 1:Nc, cf2 in 1:Nc
                
                # Skip if indices are identical (no self-interaction)
                if (ci1 == ci2 && ki1 == ki2) || (cf1 == cf2 && kf1 == kf2)
                    continue
                end
                
                # Map to global orbital indices
                i1_orb = ki1 + Nk * (ci1 - 1)
                i2_orb = ki2 + Nk * (ci2 - 1)
                f1_orb = kf1 + Nk * (cf1 - 1)
                f2_orb = kf2 + Nk * (cf2 - 1)
                
                # Apply upper triangle optimization
                if !is_upper_triangle_2body(i1_orb, i2_orb, f1_orb, f2_orb)
                    continue
                end
                
                # Calculate scattering amplitude using V_int function
                # Decompose orbital indices to get component and momentum indices
                _, ci1_comp, _ = decompose_orbital_index(i1_orb, para)
                _, ci2_comp, _ = decompose_orbital_index(i2_orb, para)
                _, cf1_comp, _ = decompose_orbital_index(f1_orb, para)
                _, cf2_comp, _ = decompose_orbital_index(f2_orb, para)
                
                # Call V_int function with proper arguments
                amp = para.V_int(kf1, kf2, ki1, ki2, cf1_comp, cf2_comp, ci1_comp, ci2_comp)
                
                # Skip zero amplitude terms
                if iszero(amp)
                    continue
                end
                
                push!(scattering_list, Scattering2Body(i1_orb, i2_orb, f1_orb, f2_orb, amp))
            end
        end
    end
    
    # Sort by (i1, i2, f1, f2) for consistent ordering
    sort!(scattering_list, by = s -> (s.i1, s.i2, s.f1, s.f2))
    
    return scattering_list
end

"""
    sort_scattering_list(scattering_list::Vector{Scattering1Body})::Vector{Scattering1Body}

Sort one-body scattering list by (i, f) indices for consistent ordering.
"""
function sort_scattering_list(scattering_list::Vector{Scattering1Body})::Vector{Scattering1Body}
    return sort(scattering_list, by = s -> (s.i, s.f))
end

"""
    sort_scattering_list(scattering_list::Vector{Scattering2Body})::Vector{Scattering2Body}

Sort two-body scattering list by (i1, i2, f1, f2) indices for consistent ordering.
"""
function sort_scattering_list(scattering_list::Vector{Scattering2Body})::Vector{Scattering2Body}
    return sort(scattering_list, by = s -> (s.i1, s.i2, s.f1, s.f2))
end

"""
    filter_zero_amplitude(scattering_list::Vector{Scattering1Body})::Vector{Scattering1Body}

Remove zero-amplitude terms from one-body scattering list.
"""
function filter_zero_amplitude(scattering_list::Vector{Scattering1Body})::Vector{Scattering1Body}
    return filter(s -> !iszero(s.amp), scattering_list)
end

"""
    filter_zero_amplitude(scattering_list::Vector{Scattering2Body})::Vector{Scattering2Body}

Remove zero-amplitude terms from two-body scattering list.
"""
function filter_zero_amplitude(scattering_list::Vector{Scattering2Body})::Vector{Scattering2Body}
    return filter(s -> !iszero(s.amp), scattering_list)
end

# Base functions for display and comparison
Base.show(io::IO, s::Scattering1Body) = print(io, "Scattering1Body: ", s.i, " → ", s.f, ", amp = ", s.amp)
Base.show(io::IO, s::Scattering2Body) = print(io, "Scattering2Body: (", s.i1, ",", s.i2, ") → (", s.f1, ",", s.f2, "), amp = ", s.amp)

Base.:(==)(a::Scattering1Body, b::Scattering1Body) = a.i == b.i && a.f == b.f && a.amp == b.amp
Base.:(==)(a::Scattering2Body, b::Scattering2Body) = a.i1 == b.i1 && a.i2 == b.i2 && a.f1 == b.f1 && a.f2 == b.f2 && a.amp == b.amp

Base.isless(a::Scattering1Body, b::Scattering1Body) = (a.i, a.f) < (b.i, b.f)
Base.isless(a::Scattering2Body, b::Scattering2Body) = (a.i1, a.i2, a.f1, a.f2) < (b.i1, b.i2, b.f1, b.f2)
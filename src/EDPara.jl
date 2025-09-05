"""
    EDPara - Parameters for momentum-conserved exact diagonalization
    
    This struct contains all parameters needed for momentum-conserved exact diagonalization
    calculations, including momentum conservation rules, component structure, and interaction
    potentials.
"""

export EDPara

"""
    @kwdef mutable struct EDPara

Parameters for momentum-conserved exact diagonalization calculations.

# Fields
- `Gk::Tuple{Int64, Int64}`: Momentum conservation mod G (default: (0, 0))
- `k_list::Matrix{Int64}`: k_list[:, i] = (k_x, k_y) momentum states
- `Nk::Int64`: Number of momentum states (default: size(k_list,2))
- `Nc_hopping::Int64`: Number of components with hopping (not conserved) (default: 1)
- `Nc_conserve::Int64`: Number of components with conserved quantum numbers (default: 1)
- `Nc::Int64`: Total number of components (default: Nc_hopping * Nc_conserve)
- `error`: Validation field that checks Nc > 0 and Nk*Nc <= 64
- `H_onebody::Array{ComplexF64,4}`: One-body Hamiltonian terms (default: zeros)
- `V_int::Function`: Interaction potential function (default: constant 1.0)

# Orbital Indexing
The orbital index formula is: i = i_k + Nk * (i_ch-1) + (Nk * Nc_hopping) * (i_cc-1)
or i = i_k + Nk * (i_c-1), where i_c = i_ch + Nc_hopping * (i_cc-1)

# Validation
- `Nc > 0`: Number of components must be positive
- `Nk*Nc <= 64`: Hilbert space dimension must not exceed 64 bits for MBS64 compatibility
"""
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
    V_int::Function
end

"""
    int_amp(i1::Int64, i2::Int64, f1::Int64, f2::Int64, para::EDPara)::ComplexF64

Rewrite V_int, inputs are (momentum + component) combined index.

Maps combined indices back to momentum and component indices, then calls the
interaction potential function with proper momentum conservation factors.
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
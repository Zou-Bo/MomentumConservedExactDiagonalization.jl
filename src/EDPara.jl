"""
    EDPara - Parameters for momentum-conserved exact diagonalization
    
    This struct contains all parameters needed for momentum-conserved exact diagonalization
    calculations, including momentum conservation rules, component structure, and interaction
    potentials.
"""

export EDPara

"""
    mutable struct EDPara

Parameters for momentum-conserved exact diagonalization calculations.

# Fields
- `Gk::Tuple{Int64, Int64}`: Momentum conservation mod G (default: (0, 0))
- `k_list::Matrix{Int64}`: k_list[:, i] = (k_x, k_y) momentum states
- `Nk::Int64`: Number of momentum states
- `Nc_hopping::Int64`: Number of components with hopping (not conserved)
- `Nc_conserve::Int64`: Number of components with conserved quantum numbers
- `Nc::Int64`: Total number of components
- `H_onebody::Array{ComplexF64,4}`: One-body Hamiltonian terms
- `V_int::Function`: Interaction potential function

# Orbital Indexing
The orbital index formula is: i = i_k + Nk * (i_ch-1) + (Nk * Nc_hopping) * (i_cc-1)
or i = i_k + Nk * (i_c-1), where i_c = i_ch + Nc_hopping * (i_cc-1)

# Validation
- `Nc > 0`: Number of components must be positive
- `Nk*Nc <= 64`: Hilbert space dimension must not exceed 64 bits for MBS64 compatibility
- `V_int` function must have correct signature: (kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) -> ComplexF64
"""
mutable struct EDPara
    # momemta are in integers

    # G integer (momentum is conserved mod G; where G=0 means no mod)
    Gk::Tuple{Int64, Int64}
    # k_list[:, i] = (k_x, k_y)
    k_list::Matrix{Int64}

    Nk::Int64
    Nc_hopping::Int64
    Nc_conserve::Int64
    Nc::Int64

    # (to be fill later by user)
    H_onebody::Array{ComplexF64,4}
    V_int::Function

    """
        EDPara(;Gk::Tuple{Int64, Int64}=(0, 0), 
               k_list::Matrix{Int64},
               Nc_hopping::Int64=1,
               Nc_conserve::Int64=1,
               H_onebody::Union{Nothing,Array{ComplexF64,4}}=nothing,
               V_int::Function)
    
    Internal constructor for EDPara with validation and default values.
    """
    function EDPara(;Gk::Tuple{Int64, Int64}=(0, 0), 
                     k_list::Matrix{Int64},
                     Nc_hopping::Int64=1,
                     Nc_conserve::Int64=1,
                     H_onebody::Union{Nothing,Array{ComplexF64,4}}=nothing,
                     V_int::Function)
        # Calculate derived fields
        Nk = size(k_list, 2)
        Nc = Nc_hopping * Nc_conserve
        
        # Validation
        @assert Nc > 0 "Number of components must be positive"
        @assert Nk*Nc <= 64 "The Hilbert space dimension must not exceed 64 bits."
        
        # Validate V_int function signature
        try
            # Test with dummy arguments to verify function signature
            test_result = V_int(1, 1, 1, 1, 1, 1, 1, 1)
            if !(test_result isa ComplexF64)
                throw(AssertionError("V_int must return ComplexF64, got $(typeof(test_result))"))
            end
        catch e
            if isa(e, AssertionError)
                rethrow(e)
            else
                throw(AssertionError("V_int function must have signature (kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) -> ComplexF64. Got error: $e"))
            end
        end
        
        # Set default H_onebody if not provided
        if H_onebody === nothing
            H_onebody = zeros(ComplexF64, Nc_hopping, Nc_hopping, Nc_conserve, Nk)
        end
        
        new(Gk, k_list, Nk, Nc_hopping, Nc_conserve, Nc, H_onebody, V_int)
    end
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
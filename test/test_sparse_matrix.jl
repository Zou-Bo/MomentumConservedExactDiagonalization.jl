using Test
using MomentumConservedExactDiagonalization
using LinearAlgebra
using SparseArrays

@testset "Sparse Matrix Construction Tests" begin
    
    @testset "Basic ED_Hmlt_Matrix functionality" begin
        # Create a simple test case
        para = EDPara(k_list=[0 1; 0 0], Nk=2, Nc_hopping=1, Nc_conserve=1)
        
        # Generate a small basis
        states = ED_mbslist(para, (1,))  # One particle
        @test length(states) == 2  # Should have 2 states for 2 orbitals, 1 particle
        
        # Create empty scattering lists for basic test
        sct1 = Vector{Scattering{1}}()
        sct2 = Vector{Scattering{2}}()
        
        # Test matrix construction
        H = ED_Hmlt_Matrix(states, sct1, sct2)
        @test size(H) == (2, 2)
        @test ishermitian(Matrix(H))
    end
    
    @testset "Threaded version consistency" begin
        # Test that threaded version produces same results
        para = EDPara(k_list=[0 1; 0 0], Nk=2, Nc_hopping=1, Nc_conserve=1)
        states = ED_mbslist(para, (1,))
        sct1 = Vector{Scattering{1}}()
        sct2 = Vector{Scattering{2}}()
        
        H_basic = ED_Hmlt_Matrix(states, sct1, sct2)
        H_threaded = ED_Hmlt_Matrix_threaded(states, sct1, sct2)
        
        # Convert to dense matrices for comparison
        @test Matrix(H_basic) ≈ Matrix(H_threaded) atol=1e-12
    end
    
    @testset "Helper functions" begin
        # Test my_searchsortedfirst
        list = [1, 3, 5, 7, 9]
        @test my_searchsortedfirst(list, 5) == 3
        @test my_searchsortedfirst(list, 4) == 0
        @test my_searchsortedfirst(list, 10) == 0
        
        # Test occ_num_between (basic test - need to create proper MBS64 state)
        # This is a placeholder as we need proper MBS64 setup
        @test_skip "occ_num_between requires proper MBS64 setup"
    end
    
    @testset "Non-empty scattering lists" begin
        # Test with actual scattering terms
        para = EDPara(k_list=[0 1; 0 0], Nk=2, Nc_hopping=1, Nc_conserve=1)
        states = ED_mbslist(para, (1,))
        
        # Create a simple one-body scattering term
        sct1 = [Scattering(1.0, 1, 2)]  # c†_1 c_2
        sct2 = Vector{Scattering{2}}()
        
        H = ED_Hmlt_Matrix(states, sct1, sct2)
        @test size(H) == (2, 2)
        @test !iszero(H)  # Should have non-zero elements
        @test ishermitian(Matrix(H))
    end
    
end
using Test
using MomentumConservedExactDiagonalization
using LinearAlgebra
using SparseArrays

@testset "MomentumConservedExactDiagonalization.jl" begin
    # Basic package loading tests
    @testset "Package loading" begin
        @test isdefined(MomentumConservedExactDiagonalization, :EDPara)
        @test isdefined(MomentumConservedExactDiagonalization, :HilbertSpace)
        @test isdefined(MomentumConservedExactDiagonalization, :HamiltonianBuilder)
        @test isdefined(MomentumConservedExactDiagonalization, :Diagonalizer)
    end

    # Test that types can be instantiated (even if empty for now)
    @testset "Type instantiation" begin
        @test_nowarn MomentumConservedExactDiagonalization.EDPara()
        @test_nowarn MomentumConservedExactDiagonalization.HilbertSpace()
        @test_nowarn MomentumConservedExactDiagonalization.HamiltonianBuilder()
        @test_nowarn MomentumConservedExactDiagonalization.Diagonalizer()
    end

    # Basic functionality tests - these will be expanded as features are implemented
    @testset "Basic functionality" begin
        # Test that we can create instances of the exported types
        para = MomentumConservedExactDiagonalization.EDPara()
        hilbert = MomentumConservedExactDiagonalization.HilbertSpace()
        hamiltonian = MomentumConservedExactDiagonalization.HamiltonianBuilder()
        diagonalizer = MomentumConservedExactDiagonalization.Diagonalizer()
        
        # Verify they are the correct types
        @test para isa MomentumConservedExactDiagonalization.EDPara
        @test hilbert isa MomentumConservedExactDiagonalization.HilbertSpace
        @test hamiltonian isa MomentumConservedExactDiagonalization.HamiltonianBuilder
        @test diagonalizer isa MomentumConservedExactDiagonalization.Diagonalizer
    end

    # Test dependencies are available
    @testset "Dependencies" begin
        # Test that LinearAlgebra and SparseArrays are available in test environment
        @test isdefined(Main, :LinearAlgebra)
        @test isdefined(Main, :SparseArrays)
        
        # Test that basic linear algebra operations work
        A = [1 2; 3 4]
        @test LinearAlgebra.det(A) ≈ -2.0
        
        # Test that sparse arrays work
        S = SparseArrays.sparse([1, 2], [1, 2], [1.0, 2.0])
        @test size(S) == (2, 2)
        @test SparseArrays.nnz(S) == 2
    end
end
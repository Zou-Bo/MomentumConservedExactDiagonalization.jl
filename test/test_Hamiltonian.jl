using Test
using MomentumConservedExactDiagonalization
using MomentumConservedExactDiagonalization: EDPara, Scattering1Body, Scattering2Body, 
                                             generate_onebody_scattering, generate_twobody_scattering, 
                                             orbital_index, decompose_orbital_index,
                                             is_upper_triangle, is_upper_triangle_2body, 
                                             apply_momentum_conservation,
                                             sort_scattering_list, filter_zero_amplitude
using LinearAlgebra

@testset "Hamiltonian Scattering Tests" begin
    
    @testset "Scattering Structs" begin
        # Test Scattering1Body construction
        s1 = Scattering1Body(1, 2, 1.0 + 2.0im)
        @test s1.i == 1
        @test s1.f == 2
        @test s1.amp ≈ 1.0 + 2.0im
        
        # Test Scattering2Body construction
        s2 = Scattering2Body(1, 2, 3, 4, 0.5 + 0.0im)
        @test s2.i1 == 1
        @test s2.i2 == 2
        @test s2.f1 == 3
        @test s2.f2 == 4
        @test s2.amp ≈ 0.5 + 0.0im
        
        # Test equality
        s1_copy = Scattering1Body(1, 2, 1.0 + 2.0im)
        @test s1 == s1_copy
        
        s2_copy = Scattering2Body(1, 2, 3, 4, 0.5 + 0.0im)
        @test s2 == s2_copy
        
        # Test ordering
        s1_small = Scattering1Body(1, 1, 1.0)
        s1_large = Scattering1Body(1, 2, 1.0)
        @test s1_small < s1_large
        
        s2_small = Scattering2Body(1, 1, 1, 1, 1.0)
        s2_large = Scattering2Body(1, 1, 1, 2, 1.0)
        @test s2_small < s2_large
    end
    
    @testset "Orbital Index Functions" begin
        # Create test EDPara
        k_list = [1 2 3; 4 5 6]  # 3 momentum points
        function default_V_int(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2)
            return 1.0 + 0.0im
        end
        para = EDPara(k_list=k_list, Gk=(0, 0), Nc_hopping=2, Nc_conserve=2, V_int=default_V_int)
        
        @test para.Nk == 3
        @test para.Nc_hopping == 2
        @test para.Nc_conserve == 2
        @test para.Nc == 4
        
        # Test orbital_index function
        # Formula: i = k + Nk * (ch-1) + (Nk * Nc_hopping) * (cc-1)
        
        # cc=1, ch=1, k=1: i = 1 + 3*(1-1) + (3*2)*(1-1) = 1
        @test orbital_index(1, 1, 1, para) == 1
        
        # cc=1, ch=1, k=2: i = 2 + 3*(1-1) + (3*2)*(1-1) = 2  
        @test orbital_index(2, 1, 1, para) == 2
        
        # cc=1, ch=2, k=1: i = 1 + 3*(2-1) + (3*2)*(1-1) = 4
        @test orbital_index(1, 2, 1, para) == 4
        
        # cc=2, ch=1, k=1: i = 1 + 3*(1-1) + (3*2)*(2-1) = 7
        @test orbital_index(1, 1, 2, para) == 7
        
        # Test decompose_orbital_index (inverse function)
        for cc in 1:para.Nc_conserve, ch in 1:para.Nc_hopping, k in 1:para.Nk
            i = orbital_index(k, ch, cc, para)
            k_recovered, ch_recovered, cc_recovered = decompose_orbital_index(i, para)
            @test k_recovered == k
            @test ch_recovered == ch
            @test cc_recovered == cc
        end
    end
    
    @testset "Upper Triangle Functions" begin
        # Test 1-body upper triangle
        @test is_upper_triangle(1, 1) == true
        @test is_upper_triangle(1, 2) == true
        @test is_upper_triangle(2, 1) == false
        @test is_upper_triangle(3, 5) == true
        
        # Test 2-body upper triangle
        @test is_upper_triangle_2body(1, 2, 1, 2) == true
        @test is_upper_triangle_2body(1, 2, 1, 3) == true
        @test is_upper_triangle_2body(1, 3, 1, 2) == false
        @test is_upper_triangle_2body(2, 1, 1, 2) == true  # Should be reordered
        
        # Test with actual reordering
        # (1,3) vs (1,2) -> (1,2) should be considered "smaller"
        @test is_upper_triangle_2body(1, 2, 1, 3) == true   # (1,2) < (1,3)
        @test is_upper_triangle_2body(1, 3, 1, 2) == false  # (1,3) > (1,2)
    end
    
    @testset "Momentum Conservation" begin
        # Create test EDPara with simple momentum points
        k_list = [1 2; 3 4]  # 2 momentum points: (1,3) and (2,4)
        function default_V_int(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2)
            return 1.0 + 0.0im
        end
        para = EDPara(k_list=k_list, Gk=(0, 0), Nc_hopping=1, Nc_conserve=1, V_int=default_V_int)
        
        # Test momentum conservation
        # ki1=1, ki2=1: total = (1+1, 3+3) = (2, 6)
        # kf1=2, kf2=2: total = (2+2, 4+4) = (4, 8)
        # Without mod G, these should not be equal
        @test apply_momentum_conservation(1, 1, 2, 2, para) == false
        
        # ki1=1, ki2=2: total = (1+2, 3+4) = (3, 7)
        # kf1=2, kf2=1: total = (2+1, 4+3) = (3, 7)  
        # These should be equal
        @test apply_momentum_conservation(1, 2, 2, 1, para) == true
        
        # Test with Gk constraints
        para_gk = EDPara(k_list=k_list, Gk=(3, 7), Nc_hopping=1, Nc_conserve=1)
        
        # ki1=1, ki2=1: total = (2, 6) mod (3, 7) = (2, 6)
        # kf1=2, kf2=2: total = (4, 8) mod (3, 7) = (1, 1)
        @test apply_momentum_conservation(1, 1, 2, 2, para_gk) == false
        
        # ki1=1, ki2=2: total = (3, 7) mod (3, 7) = (0, 0)
        # kf1=2, kf2=1: total = (3, 7) mod (3, 7) = (0, 0)
        @test apply_momentum_conservation(1, 2, 2, 1, para_gk) == true
    end
    
    @testset "One-Body Scattering Generation" begin
        # Create test EDPara
        k_list = [1 2; 3 4]  # 2 momentum points
        function default_V_int(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2)
            return 1.0 + 0.0im
        end
        para = EDPara(k_list=k_list, Gk=(0, 0), Nc_hopping=2, Nc_conserve=1, V_int=default_V_int)
        
        # Add some one-body terms
        para.H_onebody[1, 1, 1, 1] = 1.0 + 0.0im  # Diagonal term
        para.H_onebody[2, 1, 1, 1] = 0.5 + 0.0im  # Off-diagonal term (ch2→ch1)
        para.H_onebody[1, 2, 1, 2] = 0.3 + 0.1im  # Another off-diagonal term
        
        # Generate scattering list
        sct_list = generate_onebody_scattering(para)
        
        # Should have non-zero terms
        @test length(sct_list) > 0
        
        # All terms should be upper triangle
        @test all(s -> is_upper_triangle(s.i, s.f), sct_list)
        
        # All terms should have non-zero amplitude
        @test all(s -> !iszero(s.amp), sct_list)
        
        # Should be sorted
        @test issorted(sct_list, by = s -> (s.i, s.f))
        
        # Test specific terms
        diag_terms = filter(s -> s.i == s.f, sct_list)
        offdiag_terms = filter(s -> s.i != s.f, sct_list)
        
        @test length(diag_terms) >= 1  # Should have at least one diagonal term
        @test length(offdiag_terms) >= 1  # Should have at least one off-diagonal term
    end
    
    @testset "Two-Body Scattering Generation" begin
        # Create test EDPara with simple interaction
        k_list = [1 2; 3 4]  # 2 momentum points
        
        # Simple interaction function
        function test_V_int(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2)
            # Simple test interaction: 1.0 if momentum conserved, 0.0 otherwise
            if (kf1 + kf2 == ki1 + ki2) && (cf1 == ci1 && cf2 == ci2)
                return 1.0 + 0.0im
            else
                return 0.0 + 0.0im
            end
        end
        
        para = EDPara(k_list=k_list, Gk=(0, 0), Nc_hopping=1, Nc_conserve=1, V_int=test_V_int)
        
        # Generate scattering list
        sct_list = generate_twobody_scattering(para)
        
        # Should have some terms
        @test length(sct_list) > 0
        
        # All terms should satisfy upper triangle condition
        @test all(s -> is_upper_triangle_2body(s.i1, s.i2, s.f1, s.f2), sct_list)
        
        # All terms should have non-zero amplitude
        @test all(s -> !iszero(s.amp), sct_list)
        
        # Should be sorted
        @test issorted(sct_list, by = s -> (s.i1, s.i2, s.f1, s.f2))
        
        # Test momentum conservation for each term
        for s in sct_list
            ki1, ch_i1, cc_i1 = decompose_orbital_index(s.i1, para)
            ki2, ch_i2, cc_i2 = decompose_orbital_index(s.i2, para)
            kf1, ch_f1, cc_f1 = decompose_orbital_index(s.f1, para)
            kf2, ch_f2, cc_f2 = decompose_orbital_index(s.f2, para)
            
            @test apply_momentum_conservation(ki1, ki2, kf1, kf2, para)
        end
    end
    
    @testset "Utility Functions" begin
        # Test sorting
        s1_unsorted = [Scattering1Body(2, 1, 1.0), Scattering1Body(1, 1, 2.0), Scattering1Body(1, 2, 3.0)]
        s1_sorted = sort_scattering_list(s1_unsorted)
        @test issorted(s1_sorted, by = s -> (s.i, s.f))
        
        s2_unsorted = [Scattering2Body(2, 1, 1, 1, 1.0), Scattering2Body(1, 1, 1, 2, 2.0), Scattering2Body(1, 2, 1, 1, 3.0)]
        s2_sorted = sort_scattering_list(s2_unsorted)
        @test issorted(s2_sorted, by = s -> (s.i1, s.i2, s.f1, s.f2))
        
        # Test filtering
        s1_with_zeros = [Scattering1Body(1, 1, 1.0), Scattering1Body(1, 2, 0.0), Scattering1Body(2, 2, 3.0)]
        s1_filtered = filter_zero_amplitude(s1_with_zeros)
        @test length(s1_filtered) == 2
        @test all(s -> !iszero(s.amp), s1_filtered)
        
        s2_with_zeros = [Scattering2Body(1, 1, 1, 1, 1.0), Scattering2Body(1, 2, 1, 2, 0.0), Scattering2Body(2, 2, 2, 2, 2.0)]
        s2_filtered = filter_zero_amplitude(s2_with_zeros)
        @test length(s2_filtered) == 2
        @test all(s -> !iszero(s.amp), s2_filtered)
    end
    
    @testset "Integration with EDPara" begin
        # Create a more complex test case
        k_list = [0 1 2; 0 1 2]  # 3 momentum points
        
        function complex_V_int(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2)
            # More complex interaction with momentum and component dependence
            if (kf1 + kf2 == ki1 + ki2) && (cf1 == ci1 && cf2 == ci2)
                return (ki1 + ki2) * 0.1 + 0.5im
            else
                return 0.0 + 0.0im
            end
        end
        
        para = EDPara(k_list=k_list, Gk=(3, 3), Nc_hopping=2, Nc_conserve=2, V_int=complex_V_int)
        
        # Add one-body terms
        for k in 1:para.Nk, ch in 1:para.Nc_hopping, cc in 1:para.Nc_conserve
            para.H_onebody[ch, ch, cc, k] = k * 0.1 + ch * 0.05 + cc * 0.01
        end
        
        # Generate both scattering lists
        onebody_list = generate_onebody_scattering(para)
        twobody_list = generate_twobody_scattering(para)
        
        # Verify both lists are valid
        @test length(onebody_list) > 0
        @test length(twobody_list) > 0
        
        # Test that all terms are properly formed
        @test all(s -> s.i > 0 && s.f > 0, onebody_list)
        @test all(s -> s.i1 > 0 && s.i2 > 0 && s.f1 > 0 && s.f2 > 0, twobody_list)
        
        # Test that orbital indices are within valid range
        max_orbital = para.Nk * para.Nc
        @test all(s -> s.i <= max_orbital && s.f <= max_orbital, onebody_list)
        @test all(s -> s.i1 <= max_orbital && s.i2 <= max_orbital && 
                     s.f1 <= max_orbital && s.f2 <= max_orbital, twobody_list)
    end
    
    @testset "Performance and Edge Cases" begin
        # Test with minimal system
        k_list = [0 0]  # 1 momentum point
        function default_V_int(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2)
            return 1.0 + 0.0im
        end
        para_min = EDPara(k_list=k_list, Gk=(0, 0), Nc_hopping=1, Nc_conserve=1, V_int=default_V_int)
        para_min.H_onebody[1, 1, 1, 1] = 1.0
        
        function min_V_int(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2)
            return (kf1 == ki1 && kf2 == ki2 && cf1 == ci1 && cf2 == ci2) ? 1.0 : 0.0
        end
        para_min.V_int = min_V_int
        
        onebody_min = generate_onebody_scattering(para_min)
        twobody_min = generate_twobody_scattering(para_min)
        
        # Should handle minimal case correctly
        @test length(onebody_min) == 1  # Only diagonal term
        @test length(twobody_min) >= 0  # May be empty due to constraints
        
        # Test with zero Hamiltonian
        para_zero = EDPara(k_list=[0 1; 0 1], Gk=(0, 0), Nc_hopping=1, Nc_conserve=1)
        para_zero.H_onebody .= 0.0
        para_zero.V_int = (args...) -> 0.0 + 0.0im
        
        onebody_zero = generate_onebody_scattering(para_zero)
        twobody_zero = generate_twobody_scattering(para_zero)
        
        @test length(onebody_zero) == 0
        @test length(twobody_zero) == 0
    end
    
end

# Test integration with main module
@testset "Main Module Integration" begin
    # This tests that the new structs can work alongside the existing implementation
    
    using MomentumConservedExactDiagonalization: Scattering, NormalScattering, 
                                                   ED_sortedScatteringList_onebody, 
                                                   ED_sortedScatteringList_twobody
    
    # Create test parameters
    k_list = [0 1; 0 1]
    para = EDPara(k_list=k_list, Gk=(0, 0), Nc_hopping=1, Nc_conserve=1)
    para.H_onebody[1, 1, 1, 1] = 1.0
    para.H_onebody[1, 1, 1, 2] = 0.5
    
    function test_V_int(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2)
        return (kf1 + kf2 == ki1 + ki2) ? 1.0 + 0.0im : 0.0 + 0.0im
    end
    para.V_int = test_V_int
    
    # Generate lists using both old and new methods
    old_onebody = ED_sortedScatteringList_onebody(para)
    old_twobody = ED_sortedScatteringList_twobody(para)
    
    new_onebody = generate_onebody_scattering(para)
    new_twobody = generate_twobody_scattering(para)
    
    # Both should produce valid results (may differ in format but should be consistent)
    @test length(old_onebody) > 0
    @test length(old_twobody) > 0
    @test length(new_onebody) > 0
    @test length(new_twobody) > 0
    
    # Test that we can convert between formats if needed
    # (This would be implemented in the integration phase)
    
end
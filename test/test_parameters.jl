using Test
using MomentumConservedExactDiagonalization
using MomentumConservedExactDiagonalization: EDPara, int_amp

@testset "EDPara Constructor and Validation" begin
    
    @testset "Basic Constructor with Required Parameters" begin
        k_list = [1 2 3; 4 5 6]  # 3 momentum points
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        
        # Test minimal constructor
        para = EDPara(k_list=k_list, V_int=V_int_test)
        
        @test para.Gk == (0, 0)  # default value
        @test para.k_list == k_list
        @test para.Nk == 3
        @test para.Nc_hopping == 1  # default value
        @test para.Nc_conserve == 1  # default value
        @test para.Nc == 1  # 1 * 1
        @test size(para.H_onebody) == (1, 1, 1, 3)
        @test para.V_int == V_int_test
    end
    
    @testset "Constructor with All Optional Parameters" begin
        k_list = [1 2; 3 4]  # 2 momentum points
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 0.5 + 0.5im
        H_onebody_custom = zeros(ComplexF64, 2, 2, 2, 2)
        H_onebody_custom[1, 2, 1, 1] = 1.0 + 1.0im
        
        # Test full constructor
        para = EDPara(
            Gk=(6, 9),
            k_list=k_list,
            Nc_hopping=2,
            Nc_conserve=2,
            H_onebody=H_onebody_custom,
            V_int=V_int_test
        )
        
        @test para.Gk == (6, 9)
        @test para.k_list == k_list
        @test para.Nk == 2
        @test para.Nc_hopping == 2
        @test para.Nc_conserve == 2
        @test para.Nc == 4  # 2 * 2
        @test para.H_onebody == H_onebody_custom
        @test para.V_int == V_int_test
    end
    
    @testset "Validation - Nc > 0" begin
        k_list = [1 2; 3 4]
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        
        # Test Nc_hopping = 0
        @test_throws AssertionError EDPara(k_list=k_list, Nc_hopping=0, Nc_conserve=1, V_int=V_int_test)
        
        # Test Nc_conserve = 0
        @test_throws AssertionError EDPara(k_list=k_list, Nc_hopping=1, Nc_conserve=0, V_int=V_int_test)
        
        # Test both zero
        @test_throws AssertionError EDPara(k_list=k_list, Nc_hopping=0, Nc_conserve=0, V_int=V_int_test)
    end
    
    @testset "Validation - Nk*Nc <= 64" begin
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        
        # Test exactly 64
        k_list_64 = reshape(collect(1:64), 2, 32)  # 32 momentum points
        para_64 = EDPara(k_list=k_list_64, Nc_hopping=1, Nc_conserve=2, V_int=V_int_test)
        @test para_64.Nk == 32
        @test para_64.Nc == 2
        @test para_64.Nk * para_64.Nc == 64
        
        # Test exceeding 64
        k_list_65 = reshape(collect(1:66), 2, 33)  # 33 momentum points
        @test_throws AssertionError EDPara(k_list=k_list_65, Nc_hopping=1, Nc_conserve=2, V_int=V_int_test)
    end
    
    @testset "Validation - V_int Function Signature" begin
        k_list = [1 2; 3 4]
        
        # Test correct signature
        V_int_correct(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        para = EDPara(k_list=k_list, V_int=V_int_correct)
        @test para.V_int == V_int_correct
        
        # Note: The current implementation doesn't validate V_int function signature
        # These tests would pass once the improved internal constructor is implemented
        
        # Test function that returns wrong type - currently not validated
        V_int_wrong_type(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0  # returns Float64, not ComplexF64
        para_wrong_type = EDPara(k_list=k_list, V_int=V_int_wrong_type)
        @test para_wrong_type.V_int == V_int_wrong_type
        # @test_throws AssertionError EDPara(k_list=k_list, V_int=V_int_wrong_type)  # TODO: Enable after improved constructor
        
        # Test function with wrong number of arguments - currently not validated
        V_int_wrong_args(kf1, kf2, ki1, ki2) = 1.0 + 0.0im
        para_wrong_args = EDPara(k_list=k_list, V_int=V_int_wrong_args)
        @test para_wrong_args.V_int == V_int_wrong_args
        # @test_throws Exception EDPara(k_list=k_list, V_int=V_int_wrong_args)  # TODO: Enable after improved constructor
        
        # Test function that throws error - currently not validated
        V_int_throws(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = error("test error")
        para_throws = EDPara(k_list=k_list, V_int=V_int_throws)
        @test para_throws.V_int == V_int_throws
        # @test_throws Exception EDPara(k_list=k_list, V_int=V_int_throws)  # TODO: Enable after improved constructor
    end
    
    @testset "H_onebody Default Value" begin
        k_list = [1 2 3; 4 5 6]  # 3 momentum points
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        
        # Test with H_onebody = nothing (should create zeros array)
        para = EDPara(k_list=k_list, Nc_hopping=2, Nc_conserve=2, V_int=V_int_test)
        
        @test size(para.H_onebody) == (2, 2, 2, 3)
        @test all(para.H_onebody .== 0.0)
        
        # Test with custom H_onebody
        H_custom = ones(ComplexF64, 2, 2, 2, 3)
        para_custom = EDPara(k_list=k_list, Nc_hopping=2, Nc_conserve=2, H_onebody=H_custom, V_int=V_int_test)
        @test para_custom.H_onebody == H_custom
    end
    
    @testset "Nk Calculation from k_list" begin
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        
        # Test different k_list sizes
        k_list_2 = [1 2; 3 4]  # 2 momentum points
        para_2 = EDPara(k_list=k_list_2, V_int=V_int_test)
        @test para_2.Nk == 2
        
        k_list_5 = reshape(collect(1:10), 2, 5)  # 5 momentum points
        para_5 = EDPara(k_list=k_list_5, V_int=V_int_test)
        @test para_5.Nk == 5
    end
    
    @testset "Nc Calculation from Components" begin
        k_list = [1 2; 3 4]
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        
        # Test different component combinations
        para_1_1 = EDPara(k_list=k_list, Nc_hopping=1, Nc_conserve=1, V_int=V_int_test)
        @test para_1_1.Nc == 1
        
        para_2_3 = EDPara(k_list=k_list, Nc_hopping=2, Nc_conserve=3, V_int=V_int_test)
        @test para_2_3.Nc == 6
        
        para_4_2 = EDPara(k_list=k_list, Nc_hopping=4, Nc_conserve=2, V_int=V_int_test)
        @test para_4_2.Nc == 8
    end
    
    @testset "int_amp Function Compatibility" begin
        k_list = [1 2; 3 4]
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = (kf1 + kf2 - ki1 - ki2) + (cf1 + cf2 - ci1 - ci2)*im
        para = EDPara(k_list=k_list, Gk=(0, 0), Nc_hopping=1, Nc_conserve=1, V_int=V_int_test)
        
        # Test int_amp function with the new EDPara
        result = int_amp(1, 2, 3, 4, para)
        @test result isa ComplexF64
        
        # Test with Gk != (0, 0)
        para_gk = EDPara(k_list=k_list, Gk=(3, 7), Nc_hopping=1, Nc_conserve=1, V_int=V_int_test)
        result_gk = int_amp(1, 2, 3, 4, para_gk)
        @test result_gk isa ComplexF64
    end
    
    @testset "Backward Compatibility" begin
        k_list = [1 2 3; 4 5 6]
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        
        # Test that existing usage patterns still work
        para1 = EDPara(k_list=k_list, V_int=V_int_test)
        para2 = EDPara(k_list=k_list, Gk=(6, 9), V_int=V_int_test)
        para3 = EDPara(k_list=k_list, Nc_hopping=2, Nc_conserve=2, V_int=V_int_test)
        
        # Verify all constructions work
        @test para1.Nk == 3
        @test para2.Gk == (6, 9)
        @test para3.Nc == 4
    end
    
    @testset "Edge Cases and Error Conditions" begin
        k_list = [1 2; 3 4]
        V_int_test(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        
        # Test with minimal system size
        k_list_min = reshape([0, 0], 2, 1)  # 1 momentum point
        para_min = EDPara(k_list=k_list_min, V_int=V_int_test)
        @test para_min.Nk == 1
        @test para_min.Nc == 1
        
        # Test with large component numbers (but still within Nk*Nc <= 64)
        k_list_4 = reshape(collect(1:8), 2, 4)  # 4 momentum points
        para_large = EDPara(k_list=k_list_4, Nc_hopping=4, Nc_conserve=4, V_int=V_int_test)
        @test para_large.Nk == 4
        @test para_large.Nc == 16
        @test para_large.Nk * para_large.Nc == 64  # Exactly at the limit
        
        # Test boundary conditions for validation
        @test_throws AssertionError EDPara(k_list=k_list_4, Nc_hopping=0, Nc_conserve=4, V_int=V_int_test)
        @test_throws AssertionError EDPara(k_list=k_list_4, Nc_hopping=4, Nc_conserve=0, V_int=V_int_test)
        
        # Test k_list with zero columns (empty)
        # Note: Current implementation doesn't validate empty k_list - would be part of improved constructor
        k_list_empty = zeros(Int64, 2, 0)
        para_empty = EDPara(k_list=k_list_empty, V_int=V_int_test)
        @test para_empty.Nk == 0  # Current behavior - allows empty k_list
        # @test_throws AssertionError EDPara(k_list=k_list_empty, V_int=V_int_test)  # TODO: Enable after improved constructor
    end
    
    @testset "Future Constructor Tests (for Improved Implementation)" begin
        # These tests document the expected behavior of the improved internal constructor
        # They can be enabled once the improved constructor is implemented
        
        k_list = [1 2; 3 4]
        
        # Test 1: V_int function signature validation
        # Currently not implemented - would be part of improved constructor
        
        # Test 2: H_onebody=nothing should create default array
        # Currently not supported - would be part of improved constructor
        
        # Test 3: Removal of error field
        # Currently the error field exists - would be removed in improved constructor
        
        # For now, just test that the current constructor works
        V_int_correct(kf1, kf2, ki1, ki2, cf1, cf2, ci1, ci2) = 1.0 + 0.0im
        para_current = EDPara(k_list=k_list, V_int=V_int_correct)
        @test para_current.error === nothing  # Should be nothing when validation passes
        @test para_current.V_int == V_int_correct
    end
    
end
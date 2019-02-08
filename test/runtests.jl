using Correlation
using Test

using LinearAlgebra
using Statistics: mean, median
using Missings

@testset "Correlation.jl" begin
    @testset "Event Synchronisation" begin
        P = [
            0.0315    0.2412    0.2273   -0.0581    0.1386
            0.4191    0.0264   -0.3436   -0.4867    0.4199
           -0.4841    0.4578    0.3490   -0.3520    0.2408
           -0.2927    0.4860   -0.4484    0.1179   -0.1213
            0.0828    0.2984    0.1443    0.4336   -0.1386
        ]

        @testset "Detect Spikes" begin
            Z = Float64[
                 0     1     1     0     0
                 0     0     0    -1     1
                -1     1     0     0     0
                 0     1    -1     0     0
                 0     1     0     1     0
            ]

            out = Correlation.detectspikes(median, P)
            @test size(out) == size(Z)
            @test out == Z

            ex_out_p = Bool[
                 0     1     1     0     0
                 1     0     0     0     1
                 0     1     0     0     0
                 0     1     0     0     0
                 0     1     0     1     0
            ]
            ex_out_n = Bool[
                 0     0     0     0     0
                 0     0     0     1     0
                 1     0     0     0     0
                 1     0     1     0     0
                 0     0     0     0     0
            ]

            out = Correlation.threshold_p(mean, P)
            @test ex_out_p == out
            out = Correlation.threshold_n(mean, P)
            @test ex_out_n == out
            @test_throws MethodError Correlation.detectspikes(identity, P)
        end

        @testset "escor" begin
            out = Correlation.sqrtescor(median, P, 1)
            expected_out = [
                0.8201   -0.1917    0.1917    0.0000    0.5038
               -0.1917    0.6848   -0.6848    0.0000   -0.1591
                0.1917   -0.6848    0.6848   -0.0000    0.1591
                0.0000    0.0000    0.0000    1.0000    0.0000
                0.5038   -0.1591    0.1591    0.0000    0.8340
            ]
            @test size(out) == size(expected_out)
            @test isapprox(out, expected_out; rtol=1e-4)

            @testset "zero on diagonal" begin
                # https://gitlab.invenia.ca/invenia/Simulation.jl/merge_requests/83#note_47700
                @test Correlation.sqrtescor(median, Diagonal([1, 0]), 1) ≈ [1 0; 0 0]
            end
        end

        @testset "Missing data" begin
            # To test against MATLAB: out_mean = escorr(P, 1, 'mean')
            # where P is the following with missing -> NaN

            P = [
                0.1443   -0.2923   -0.1889    0.0949   -0.4145
               -0.1214   missing    0.4234   -0.2378   -0.2375
                0.3116   -0.0291   -0.0698    0.1028    0.3010
                0.0328   -0.2695   -0.3152    0.2112   -0.4708
               -0.1493    0.3443    0.4049   -0.2783    0.4289
                0.4390   -0.3052    0.4797   -0.3826    0.2303
                0.3759   -0.2741   -0.0611   -0.2033   -0.0114
                0.0502   -0.3293   -0.3889   -0.1812    0.0785
                0.1225   -0.2723   -0.2419   missing   -0.2627
                0.0870   -0.0643   -0.0913    0.0079   -0.0412
            ]

            out_mean = Correlation.sqrtescor(mean, P, 1)
            out_median = Correlation.sqrtescor(median, P, 1)

            expected_out_mean = [
                0.6536    0.3456   -0.0377    0.0605   -0.0045    0.4187   -0.0221   -0.3929   -0.0221    0.3429
                0.3456    0.5483   -0.0557   -0.3048   -0.1396    0.1934    0.4040   -0.0151    0.4040    0.3167
               -0.0377   -0.0557    0.7246   -0.3630   -0.0344   -0.0035   -0.2148    0.1846   -0.2148    0.4594
                0.0605   -0.3048   -0.3630    0.6917   -0.1905   -0.0580   -0.2170   -0.2956   -0.2170   -0.2685
               -0.0045   -0.1396   -0.0344   -0.1905    0.7818    0.2459   -0.2974   -0.1177   -0.2974   -0.2840
                0.4187    0.1934   -0.0035   -0.0580    0.2459    0.6369   -0.1009   -0.5225   -0.1009    0.1563
               -0.0221    0.4040   -0.2148   -0.2170   -0.2974   -0.1009    0.5486    0.1959    0.5486    0.0639
               -0.3929   -0.0151    0.1846   -0.2956   -0.1177   -0.5225    0.1959    0.6003    0.1959   -0.0012
               -0.0221    0.4040   -0.2148   -0.2170   -0.2974   -0.1009    0.5486    0.1959    0.5486    0.0639
                0.3429    0.3167    0.4594   -0.2685   -0.2840    0.1563    0.0639   -0.0012    0.0639    0.6211
            ]

            expected_out_median = [
                0.5674    0.3160    0.3300    0.0410   -0.0410    0.3692    0.0004   -0.3252   -0.3392    0.3300
                0.3160    0.5560    0.2906   -0.2095    0.2095    0.2391    0.4850   -0.0743    0.1911    0.2906
                0.3300    0.2906    0.5735   -0.2162    0.2162   -0.0828    0.0584    0.1924   -0.0904    0.5735
                0.0410   -0.2095   -0.2162    0.5945   -0.5945   -0.0975   -0.0954   -0.2636   -0.2569   -0.2162
               -0.0410    0.2095    0.2162   -0.5945    0.5945    0.0975    0.0954    0.2636    0.2569    0.2162
                0.3692    0.2391   -0.0828   -0.0975    0.0975    0.6925   -0.0035   -0.5091   -0.1872   -0.0828
                0.0004    0.4850    0.0584   -0.0954    0.0954   -0.0035    0.7083    0.0579    0.4845    0.0584
               -0.3252   -0.0743    0.1924   -0.2636    0.2636   -0.5091    0.0579    0.5679    0.3012    0.1924
               -0.3392    0.1911   -0.0904   -0.2569    0.2569   -0.1872    0.4845    0.3012    0.5827   -0.0904
                0.3300    0.2906    0.5735   -0.2162    0.2162   -0.0828    0.0584    0.1924   -0.0904    0.5735
            ]

            @test size(out_mean) == size(expected_out_mean)
            @test isapprox(out_mean, expected_out_mean; rtol=1e-4)
            @test size(out_median) == size(expected_out_median)
            @test isapprox(out_median, expected_out_median; rtol=1e-4)
        end

        @testset "missing_columns" begin
            r = allowmissing(ones(5,4))
            @test Correlation.missing_columns(r) == falses(4)
            r[2,2] = missing
            @test Correlation.missing_columns(r) == falses(4)
            r[:,2] .= missing
            @test Correlation.missing_columns(r) == [false, true, false, false]
        end
    end
end

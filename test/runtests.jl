include("strainStressTest.jl")
include("stiffnessTensorTest.jl")
include("typetest.jl")
include("Elasticity1DLaminateTest.jl")
include("ElasticityHashinEllipsoidGenerationTest.jl")
include("ElasticityHashinEllipsoidTest.jl")





#= @testset "Misc" begin =#
#=   fromNormal = [1.2,2.3,5.1] =#
#=   toNormalA = [1.2,2.3,5.1+1e-13] =#
#=   toNormalB = [2.2,4.4,5.3] =#

#=   RA = RotationMatrixNormal(fromNormal,toNormalA) =#
#=   RB = RotationMatrixNormal(fromNormal,toNormalB) =#

#=   @test isapprox(RA'*RA,eye(3),rtol=0.0,atol=1e-12) =#
#=   @test isapprox(RB'*RB,eye(3),rtol=0.0,atol=1e-12) =#
#=   @test isapprox(det(RA),1,rtol=0.0,atol=1e-12) =#
#=   @test isapprox(det(RB),1,rtol=0.0,atol=1e-12) =#

#=   @test_throws ArgumentError DepolarizationFactors(0.0,0.0,0.0) =#
#=   @test_throws ArgumentError DepolarizationFactors(1.0,2.0,1.0) =#

#=   (d1,d2,d3) = DepolarizationFactors(2.0,2.0,2.0) =#
#=   @test isapprox(d1+d2+d3,1.0,rtol=0.0,atol=1e-12) =#

#=   (d1,d2,d3) = DepolarizationFactors(2.0,1.0,1.0) =#
#=   @test isapprox(d1+d2+d3,1.0,rtol=0.0,atol=1e-12) =#

#=   (d1,d2,d3) = DepolarizationFactors(0.5,1.0,1.0) =#
#=   @test isapprox(d1+d2+d3,1.0,rtol=0.0,atol=1e-12) =#

#=   (d1,d2,d3) = DepolarizationFactors(0.5,1.0,Inf) =#
#=   @test isapprox(d1+d2+d3,1.0,rtol=0.0,atol=1e-12) =#
#= end =#

#= #@testset "Geometries" begin =#
#= #  M = Lattice([8;8;1]) =#
#= #  (C,epsilon0) = Elasticity_Hashin_ellipsoid(M) =#
#= #  (epsilon_analytic,tau0) = Elasticity_Hashin_ellipsoid_analytic(M) =#
#= # =#
#= #end =#

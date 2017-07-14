include("typetest.jl")
include("Elasticity1DLaminateTest.jl")
include("ElasticityHashinEllipsoidGenerationTest.jl")
include("ElasticityHashinEllipsoidTest.jl")


#= @testset "Strain and Stress" begin =#
#=   strain = Strain(Float64) =#
#=   for i in 1:6 =#
#=     strain.val[i] = 1.0i =#
#=     @test strain[i] == strain.val[i] =#
#=   end =#
#=   @test strain[1,1] == strain.val[1] =#
#=   @test 2.0strain[1,2] == strain.val[6] =#
#=   @test 2.0strain[1,3] == strain.val[5] =#
#=   @test 2.0strain[2,1] == strain.val[6] =#
#=   @test strain[2,2] == strain.val[2] =#
#=   @test 2.0strain[2,3] == strain.val[4] =#
#=   @test 2.0strain[3,1] == strain.val[5] =#
#=   @test 2.0strain[3,2] == strain.val[4] =#
#=   @test strain[3,3] == strain.val[3] =#

#=   for i in 1:6 =#
#=     strain[i] = 1.0i+10.0 =#
#=     @test strain.val[i] == 1.0i+10.0 =#
#=   end =#

#=   for i in 1:3 =#
#=     for j in 1:3 =#
#=       r = rand() =#
#=       strain[i,j] = r =#
#=       i == j && @test strain.val[i] == r =#
#=       i == 1 && j == 2 && @test strain.val[6] == 2.0r =#
#=       i == 1 && j == 3 && @test strain.val[5] == 2.0r =#
#=       i == 2 && j == 1 && @test strain.val[6] == 2.0r =#
#=       i == 2 && j == 3 && @test strain.val[4] == 2.0r =#
#=       i == 3 && j == 1 && @test strain.val[5] == 2.0r =#
#=       i == 3 && j == 2 && @test strain.val[4] == 2.0r =#
#=     end =#
#=   end =#


#=   stress = Stress(Float64) =#
#=   for i in 1:6 =#
#=     stress.val[i] = 1.0i =#
#=     @test stress[i] == stress.val[i] =#
#=   end =#
#=   @test stress[1,1] == stress.val[1] =#
#=   @test stress[1,2] == stress.val[6] =#
#=   @test stress[1,3] == stress.val[5] =#
#=   @test stress[2,1] == stress.val[6] =#
#=   @test stress[2,2] == stress.val[2] =#
#=   @test stress[2,3] == stress.val[4] =#
#=   @test stress[3,1] == stress.val[5] =#
#=   @test stress[3,2] == stress.val[4] =#
#=   @test stress[3,3] == stress.val[3] =#

#=   for i in 1:6 =#
#=     stress[i] = 1.0i+10.0 =#
#=     @test stress.val[i] == 1.0i+10.0 =#
#=   end =#

#=   for i in 1:3 =#
#=     for j in 1:3 =#
#=       r = rand() =#
#=       stress[i,j] = r =#
#=       i == j && @test stress.val[i] == r =#
#=       i == 1 && j == 2 && @test stress.val[6] == r =#
#=       i == 1 && j == 3 && @test stress.val[5] == r =#
#=       i == 2 && j == 1 && @test stress.val[6] == r =#
#=       i == 2 && j == 3 && @test stress.val[4] == r =#
#=       i == 3 && j == 1 && @test stress.val[5] == r =#
#=       i == 3 && j == 2 && @test stress.val[4] == r =#
#=     end =#
#=   end =#
#= end =#

#= @testset "Multiplication with stiffness tensors" begin =#
#=   stress = Stress(Float64) =#
#=   strain = Strain(Float64) =#

#=   lambda = LamesFirstParameter(0.45) =#
#=   mu = ShearModulus(20.2) =#
#=   CIsotropic = IsotropicStiffnessTensor(lambda,mu) =#

#=   strain.val = [ 9.571669482429456 ,4.853756487228412 ,8.002804688888002 =#
#=                 ,1.418863386272153 ,4.217612826262750 ,9.157355251890671] =#

#=   mult!(stress,CIsotropic,strain) =#

#=   #Computed by Matlab implementation =#
#=   @test stress.val ≈ 1.0e+02 * =#
#=   [ =#
#=    3.967881508864957, =#
#=    2.061844658803735, =#
#=    3.334060132274209, =#
#=    0.286610404026975, =#
#=    0.851957790905075, =#
#=    1.849785760881915 =#
#=   ] =#

#=   E_p = YoungsModulus(123.4) =#
#=   E_t = YoungsModulus(1001.321) =#
#=   nu_p = PoissonsRatio(0.2) =#
#=   nu_pt = PoissonsRatio(0.12) =#
#=   nu_tp = PoissonsRatio(0.32) =#
#=   mu_t = ShearModulus(22.43) =#

#=   #CTransversal = TransversalIsotropicZStiffnessTensor(E_p,E_t,nu_p,nu_pt,nu_tp,mu_t) =#

#=   #mult!(stress,CTransversal,strain) =#

#=   ##Computed by Matlab implementation =#
#=   #@test stress.val ≈ 1.0e+03 * =#
#=   #[ =#
#=   # 1.910255706466427, =#
#=   # 1.425096986793253, 9.652008207681787, =#
#=   # 0.015912552877042, =#
#=   # 0.047300527846537, =#
#=   # 0.235420341267356 =#
#=   #] =#
#= end =#

#= @testset "Gamma0 operator" begin =#
#=   stress = Stress(Complex128) =#
#=   strain = Strain(Complex128) =#
#=   gamma = Gamma0() =#

#=   patternPoints = [0.678735154857773, 0.757740130578333, 0.743132468124916] =#

#=   lambda = LamesFirstParameter(0.45) =#
#=   mu = ShearModulus(20.2) =#
#=   CIsotropic = IsotropicStiffnessTensor(lambda,mu) =#

#=   stress.val = complex([ 9.571669482429456 ,4.853756487228412 ,8.002804688888002 =#
#=                         ,1.418863386272153 ,4.217612826262750 ,9.157355251890671]) =#

#=   applyGammaHat!(strain,gamma,stress,CIsotropic,patternPoints) =#

#=   @test strain.val ≈ =# 
#=   [ =#
#=    -0.227739399421220, =#
#=    -0.105376068690819, =#
#=    -0.081609049045277, =#
#=    -0.186557864007233, =#
#=    -0.323884040352719, =#
#=    -0.348637477862576 =#
#=   ]+0.0im =#

#=   applyGammaHat!(strain,gamma,stress,CIsotropic,[0.0,0.0,0.0]) =#

#=   @test norm(strain.val) == 0.0 =#
#= end =#

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

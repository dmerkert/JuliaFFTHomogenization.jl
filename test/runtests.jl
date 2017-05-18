using JuliaBasicScheme
using Base.Test

@testset "Multiplication with stiffness tensors" begin
  stress = Stress(Float64)
  strain = Strain(Float64)

  lambda = LamesFirstParameter(0.45)
  mu = ShearModulus(20.2)
  CIsotropic = IsotropicStiffnessTensor(lambda,mu)

  strain.val = [ 9.571669482429456 ,4.853756487228412 ,8.002804688888002
                ,1.418863386272153 ,4.217612826262750 ,9.157355251890671]

  mult!(stress,CIsotropic,strain)

  #Computed by Matlab implementation
  @test isapprox(
                 stress.val,
                 1.0e+02 *
                 [
                  3.967881508864957,
                  2.061844658803735,
                  3.334060132274209,
                  0.286610404026975,
                  0.851957790905075,
                  1.849785760881915
                 ]
                )

  E_p = YoungsModulus(123.4)
  E_t = YoungsModulus(1001.321)
  nu_p = PoissonsRatio(0.2)
  nu_pt = PoissonsRatio(0.12)
  nu_tp = PoissonsRatio(0.32)
  mu_t = ShearModulus(22.43)

  CTransversal = TransversalIsotropicZ(E_p,E_t,nu_p,nu_pt,nu_tp,mu_t)

  mult!(stress,CTransversal,strain)

  #Computed by Matlab implementation
  @test isapprox(
                 stress.val,
                 1.0e+03 *
                 [
                  1.910255706466427,
                  1.425096986793253,
                  9.652008207681787,
                  0.015912552877042,
                  0.047300527846537,
                  0.235420341267356
                 ]
                )
end

@testset "Gamma0 operator" begin
  stress = Stress(Complex64)
  strain = Strain(Complex64)

  patternPoints = [0.678735154857773, 0.757740130578333, 0.743132468124916]

  lambda = LamesFirstParameter(0.45)
  mu = ShearModulus(20.2)
  CIsotropic = IsotropicStiffnessTensor(lambda,mu)

  stress.val = complex([ 9.571669482429456 ,4.853756487228412 ,8.002804688888002
                        ,1.418863386272153 ,4.217612826262750 ,9.157355251890671])

  mult!(strain,stress,CIsotropic,patternPoints)
  @show real(strain.val)
  @show [
         -0.227739399421220,
         -0.105376068690819,
         -0.081609049045277,
         -0.186557864007233,
         -0.257444439344756,
         -0.318037320319807
        ]

  #@test isapprox(strain.val,
  #               [
  #                -0.227739399421220,
  #                -0.105376068690819,
  #                -0.081609049045277,
  #                -0.186557864007233,
  #                -0.257444439344756,
  #                -0.318037320319807
  #               ]
  #              )
  
end

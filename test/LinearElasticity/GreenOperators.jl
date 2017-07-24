using JuliaFFTHomogenization
using Base.Test

@testset "Gamma0 operator" begin
  stress = Stress{Complex128}()
  strain = Strain{Complex128}()
  gamma = Gamma0()

  patternPoints = [0.678735154857773, 0.757740130578333, 0.743132468124916]

  lambda = LamesFirstParameter(0.45)
  mu = ShearModulus(20.2)
  CIsotropic = IsotropicStiffnessTensor(lambda,mu)

  stress.val = complex([ 9.571669482429456 ,4.853756487228412 ,8.002804688888002
                        ,1.418863386272153 ,4.217612826262750 ,9.157355251890671])

  mult!(strain,gamma,stress,CIsotropic,patternPoints,zeros((6,6)))

  @test strain.val ≈ 
  [
   -0.227739399421220,
   -0.105376068690819,
   -0.081609049045277,
   -0.186557864007233,
   -0.323884040352719,
   -0.348637477862576
  ]+0.0im

  mult!(strain,gamma,stress,CIsotropic,[0.0,0.0,0.0],zeros((6,6)))

  @test norm(strain.val) == 0.0
end

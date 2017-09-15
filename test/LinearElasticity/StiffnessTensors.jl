using JuliaFFTHomogenization
using Base.Test

@testset "StiffnessTensors" begin
  v = [1.0;2.0;3.0;4.0;5.0;6.0]
  strain = Strain(v)
  stress = Stress(v)

  l = LamesFirstParameter()
  mu = ShearModulus()
  CIso = IsotropicStiffnessTensor(l,mu)

  E_p = YoungsModulus()
  E_t = YoungsModulus()
  nu_p = PoissonsRatio()
  nu_pt = PoissonsRatio()
  nu_tp = PoissonsRatio()
  mu_t = ShearModulus()
  CTrans = TransversalIsotropicZStiffnessTensor(E_p,E_t,nu_p,nu_pt,nu_tp,mu_t)

  CDiag = DiagonalStiffnessTensor(v)

  CAnis = AnisotropicStiffnessTensor(diagm(v))

  JuliaFFTHomogenization.eig(CIso)
  JuliaFFTHomogenization.eig(CTrans)
  JuliaFFTHomogenization.eig(CDiag)
  JuliaFFTHomogenization.eig(CAnis)

  mult!(stress,CIso,strain)
  mult!(stress,CTrans,strain)
  mult!(stress,CDiag,strain)
  mult!(stress,CAnis,strain)

  convert(AnisotropicStiffnessTensor,CIso)
  convert(AnisotropicStiffnessTensor,CTrans)
  convert(AnisotropicStiffnessTensor,CDiag)

  convert(TransversalIsotropicZStiffnessTensor,CIso)

  for t1 in (CIso,CTrans,CDiag,CAnis)
    for t2 in (CIso,CTrans,CDiag,CAnis)
      t1+t2
      t1-t2
    end
  end
end

@testset "Composite StiffnessTensors" begin
  E1 = YoungsModulus(1.0)
  E2 = YoungsModulus(2.0)
  l1 = LamesFirstParameter(0.1)
  l2 = LamesFirstParameter(0.2)

  C1 = IsotropicStiffnessTensor(E1,l1)
  C2 = IsotropicStiffnessTensor(E2,l2)

  volumes = [0.25;0.75]

  CArithmetic = CompositeArithmeticMeanStiffnessTensor([C1;C2],volumes)
  CArithmeticA = convert(AnisotropicStiffnessTensor,CArithmetic)

  @test volumes[1]*convert(AnisotropicStiffnessTensor,C1).C + 
  volumes[2]*convert(AnisotropicStiffnessTensor,C2).C ≈
  convert(AnisotropicStiffnessTensor,CArithmeticA).C

  @test inv(C1).C*convert(AnisotropicStiffnessTensor,C1).C ≈ eye(6)

  CHarmonic = CompositeHarmonicMeanStiffnessTensor([C1;C2],volumes)
  CHarmonicA = convert(AnisotropicStiffnessTensor,CHarmonic)

  @test inv(
            volumes[1]*inv(convert(AnisotropicStiffnessTensor,C1).C) + 
            volumes[2]*inv(convert(AnisotropicStiffnessTensor,C2).C)
           ) ≈ CHarmonicA.C
  
  CArithmeticHarmonic = CompositeAvgArithmeticHarmonicMeanTensor([C1;C2],volumes)
  CArithmeticHarmonicA = convert(AnisotropicStiffnessTensor,CArithmeticHarmonic)

  @test 0.5*(CArithmeticA.C+CHarmonicA.C) ≈ CArithmeticHarmonicA.C




end

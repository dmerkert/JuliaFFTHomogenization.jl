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

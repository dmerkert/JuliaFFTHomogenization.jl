using JuliaFFTHomogenization
using Base.Test



@testset "Coefficient Tensor Fields" begin
  l = LamesFirstParameter(0.243443)
  mu = ShearModulus(6.34331247)

  CIso = IsotropicStiffnessTensor(l,mu)
  CIso2 = IsotropicStiffnessTensor(l+1,mu-1)
  CTrans = convert(TransversalIsotropicZStiffnessTensor,CIso)
  CTrans2 = convert(TransversalIsotropicZStiffnessTensor,CIso2)

  CIsoField = CoefficientTensorField{IsotropicStiffnessTensor}((64,64))
  CTransField =
  CoefficientTensorField{TransversalIsotropicZStiffnessTensor}((64,64))
  for coord in CartesianRange(size(CIsoField))
    l = LamesFirstParameter(0.5rand())
    mu = ShearModulus(rand()+1)

    CIsoField[coord] = IsotropicStiffnessTensor(l,mu)
    CTransField[coord] =
    convert(TransversalIsotropicZStiffnessTensor,CIsoField[coord])
  end
  solver = BasicScheme{CauchyConvergenceCriterion}()
  C0 = getReferenceTensor(CIsoField,solver)
  C0Trans = getReferenceTensor(CTransField,solver)

  C0An = convert(AnisotropicStiffnessTensor,C0)
  C0An2 = convert(AnisotropicStiffnessTensor,C0Trans)
end

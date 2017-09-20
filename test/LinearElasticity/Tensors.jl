using JuliaFFTHomogenization
using Base.Test

function f(i,j,k,l)
  2.0i^2 +
  2.3j^2 +
  2.6k^2 +
  2.7l^2 +
  (i*j-3.0k*l)*i^2
end

function f1(i,j,k,l)
  1.0
end

@testset "Tensors" begin
  T = Array{Float64}((6,6))
  TI = Array{Float64}((6,6))

  toVoigtSSymFourthOrder!(T,f)
  toInverseVoigtSSymFourthOrder!(TI,f)

  @test issymmetric(T)
  @test issymmetric(TI)

  toVoigtSSymFourthOrder!(T,f1)
  toInverseVoigtSSymFourthOrder!(TI,f1)

  @test all(T .== 1.0)
  @test all(TI[1:3,1:3] .== 1.0)
  @test all(TI[4:6,4:6] .== 4.0)
  @test all(TI[1:3,4:6] .== 2.0)
  @test all(TI[4:6,1:3] .== 2.0)
end

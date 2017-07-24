using JuliaFFTHomogenization
using Base.Test
using MPAWL

@testset "Transformation" begin
  transform = TruncatedTrigonometricPolynomials()

  L = Lattice(diagm([3;4;5]))

  strain = SolutionTensorField{Float64,Strain}(rand((6,L.size...)))
  strainC = SolutionTensorField{Complex128,Strain}(rand(Complex128,(6,L.size...)))

  avg = sum(strain.val,(2))[:]/(3*4*5)

  transform!(
             strainC,
             strain,
             transform,
             L
            )

  @test strainC.val[:,1,1,1]/(3*4*5) ≈ avg

  setAveragingFrequency!(
                         strainC,
                         Strain([1.0;2.0;3.0;4.0;5.0;6.0]),
                         transform,
                         L
                        )

  transformInverse!(
                    strain,
                    strainC,
                    transform,
                    L
                   )

  avg = sum(strain.val,(2))[:]/(3*4*5)

  @test [1.0;2.0;3.0;4.0;5.0;6.0] ≈ avg

end

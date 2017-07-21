using JuliaFFTHomogenization
using Base.Test

@testset "Basic Scheme" begin
  BasicScheme(
              tol = 1e-8,
              maxIter = 100,
              verbose = true,
              printSkip = 1,
              convergenceCriterion = NormConvergenceCriterion()
             )

end


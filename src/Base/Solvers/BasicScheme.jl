export BasicScheme

immutable BasicScheme <: Solver
  tol       :: Float64
  maxIter   :: Int
  verbose   :: Bool
  printSkip :: Int
  convergenceCriterion :: ConvergenceCriterion

  BasicScheme(;
                 tol :: Float64 = 1e-6,
                 maxIter :: Int = 100,
                 verbose :: Bool = false,
                 printSkip :: Int = 50,
                 convergenceCriterion :: C = CauchyConvergenceCriterion()
                ) where {C <: ConvergenceCriterion} = new(
                                  tol,
                                  maxIter,
                                  verbose,
                                  printSkip,
                                  convergenceCriterion
                                 )
end

function _solve!(
                 gradient            :: SolutionTensorField{R,G,N},
                 flux                :: SolutionTensorField{R,F,N},
                 gradientFourier     :: SolutionTensorField{C,G,N},
                 fluxFourier         :: SolutionTensorField{C,F,N},
                 coefficientField    :: CoefficientTensorField{T,M},
                 macroscopicGradient :: G,
                 solver              :: BasicScheme,
                 gamma               :: GreenOperator,
                 ansatzSpace         :: AnsatzSpace,
                 lattice             :: Lattice
                ) where {
                         G <: GradientSolutionTensor,
                         F <: FluxSolutionTensor,
                         R <: Real,
                         C <: Complex,
                         N,
                         M,
                         T
                        }

  init!(gradient,macroscopicGradient)
  referenceTensor = getReferenceTensor(coefficientField,solver)


  init!(solver.convergenceCriterion,gradient)

  error = Inf
  iterationStep = 0
  timeStart = time()
  elapedTime = 0

  if solver.verbose
    @printf "%-13s%-15s%-17s\n" "Iteration" "Distance" "Elapsed (seconds)"
    println(repeat("-", 45))
  end

  values2Coefficients!(gradient,ansatzSpace,lattice)
  while (error > solver.tol) && (iterationStep <= solver.maxIter)
   flux = mult!(flux,
                 coefficientField,
                 referenceTensor,
                 gradient
                )
    transform!(fluxFourier,
               flux,
               ansatzSpace,
               lattice
              )
    mult!(gradientFourier,
          gamma,
          fluxFourier,
          referenceTensor,
          lattice,
          ansatzSpace
         )
    setAveragingFrequency!(gradientFourier,
                           macroscopicGradient,
                           ansatzSpace,
                           lattice
                          )
    transformInverse!(gradient,
                      gradientFourier,
                      ansatzSpace,
                      lattice
                     )

    error = computeError(solver.convergenceCriterion,gradient)
    iterationStep += 1
    elapsedTime = time() - timeStart
    if solver.verbose
      @printf "%-13i%-15.5e%-18.5f\n" iterationStep error elapsedTime
    end
    set!(solver.convergenceCriterion,gradient)
  end
  coefficients2Values!(gradient,ansatzSpace,lattice)
  gradient
end

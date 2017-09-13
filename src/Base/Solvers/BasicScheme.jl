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
  elapsedTime = 0.0

  if solver.verbose
    sizeGradient = Base.summarysize(gradient)
    sizeGradientFourier = Base.summarysize(gradientFourier)
    sizeCoefficients = Base.summarysize(coefficientField)
    sizeTotal = sizeGradient + sizeGradientFourier + sizeCoefficients
    println(
            string(
                   "Size Gradient: ",round(sizeGradient/(1024^2),1),"MB\n",
                   "Size Gradient Fourier: ",round(sizeGradientFourier/(1024^2),1),"MB\n",
                   "Size Coefficients: ",round(sizeCoefficients/(1024^2),1),"MB\n",
                   "Size Total: ", round(sizeTotal/(1024^2),1),"MB\n"
                  )
           )
    @printf "%-13s%-15s%-17s\n" "Iteration" "Distance" "Elapsed (seconds)"
    println(repeat("-", 45))
  end

  values2Coefficients!(gradient,ansatzSpace,lattice)
  while (error > solver.tol) && (iterationStep <= solver.maxIter)
    #= asd = copy(gradient) =#

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
  converged = (error < solver.tol)
  (converged,elapsedTime,iterationStep,error,solver.maxIter,solver.tol)
end

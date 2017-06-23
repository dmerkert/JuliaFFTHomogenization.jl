immutable BasicScheme <: Solver
  tol       :: Float64
  maxIter   :: Int
  verbose   :: Bool
  printSkip :: Int

  BasicScheme(;
               tol=1e-6,
               maxIter=100,
               verbose=false,
               printSkip=50
              ) = new(tol,maxIter,verbose,printSkip)
end

function solve!(problem :: EffectiveTensorProblem,
                approximationMethod :: ApproximationMethod,
                solver :: Solver)

  initializeProblem!(problem)
  for i = 1:length(problem)
    subproblem = problem[i]
    solve!(subproblem,approximationMethod,solver)
    problem[i] = subproblem
  end
  problem
end

function solve!(problem :: MacroscopicGradientProblem,
                approximationMethod :: ApproximationMethod,
                solver :: Solver)
  initializeProblem!(problem)
  _solve!(unpackProblem(problem)...,
          solver,
          approximationMethod.gamma,
          approximationMethod.transformation,
          problem.lattice
         )
  postprocess!(problem)
  problem
end

function _solve!{G <: GradientSolutionTensor,
                 F <: FluxSolutionTensor,
                 R <: Real,
                 C <: Complex
                }(
                  gradient            :: SolutionTensorField{R,G},
                  flux                :: SolutionTensorField{R,F},
                  gradientFourier     :: SolutionTensorField{C,G},
                  fluxFourier         :: SolutionTensorField{C,F},
                  coefficientField    :: CoefficientTensorField,
                  macroscopicGradient :: G,
                  solver              :: BasicScheme,
                  gamma               :: GreenOperator,
                  transformation      :: Transformation,
                  lattice             :: Lattice
                   )

  init!(gradient,macroscopicGradient)
  referenceTensor = getReferenceTensor(coefficientField,solver)



  gradientPrev = copy(gradient)
  error = Inf
  iterationStep = 0
  timeStart = time()
  elapedTime = 0

  @printf "%-13s%-15s%-17s\n" "Iteration" "Distance" "Elapsed (seconds)"
  println(repeat("-", 45))

  while (error > solver.tol) && (iterationStep <= solver.maxIter)
    flux = mult!(flux,
                 coefficientField-referenceTensor,
                 gradient
                )
    transform!(fluxFourier,
               flux,
               transformation,
               lattice
              )
    mult!(gradientFourier,
          gamma,
          fluxFourier,
          referenceTensor,
          lattice
         )
    setAveragingFrequency!(gradientFourier,
                           macroscopicGradient,
                           transformation,
                           lattice
                          )
    transformInverse!(gradient,
                      gradientFourier,
                      transformation,
                      lattice
                     )

    error = norm(gradient-gradientPrev)/norm(gradient)
    iterationStep += 1
    elapsedTime = time() - timeStart
    @printf "%-13i%-15.5e%-18.5f\n" iterationStep error elapsedTime
    copy!(gradientPrev,gradient)
  end
  gradient
end

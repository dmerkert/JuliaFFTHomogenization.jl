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

  mgr = DefaultManager(solver.tol,
                       solver.maxIter,
                       solver.verbose,
                       solver.printSkip)

  istate = DefaultState(gradient)

  iii = managed_iteration(mgr, istate; by=(x,y)->norm(x-y)) do _gradient
    __gradient = copy(_gradient)
    flux = mult!(flux,
                 coefficientField-referenceTensor,
                 __gradient
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
    transformInverse!(__gradient,
                      gradientFourier,
                      transformation,
                      lattice
                     )
    __gradient
  end
  @show iii.n
  @show iii.change

end

immutable BasicScheme <: Solver
  tol :: Float64
  maxIter :: Int
  verbose :: Bool
  printSkip :: Int

  BasicScheme(;
               tol=1e-6,
               maxIter=100,
               verbose=false,
               printSkip=50
              ) = new(lattice,gamma,tol,maxIter,verbose,printSkip)
end

function solve!(problem :: EffectiveTensorProblem,
                approximationMethod :: ApproximationMethod,
                solver :: Solver)

  initializeProblem!(problem)
  for i = 1:size(problem)
    subproblem = getMacroscopicGradientProblem(problem,i)
    solve!(subproblem,approximationMethod,solver)
    setMacroscopicGradientProblem!(problem,subproblem,i)
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

function _solve!{G <: GradientSolutionTensor, R}(gradient ::
                                                 SolutionTensorField{R,G},
                coefficientField :: CoefficientTensorField,
                macroscopicGradient :: G,
                solver :: BasicScheme,
                gamma :: GreenOperator,
                transformation :: Transformation,
                lattice :: Lattice
               )

  initialize!(gradientField,macroscopicGradient)

  referenceCoefficient = zeros(33)

  mgr = DefaultManager(solver.tol,
                       solver.maxIter,
                       solver.verbose,
                       solver.printSkip)
  istate = DefaultState(gradientField)

  managed_iteration(mgr, istate; by=(x,y)->norm(x-y)) do gradient
   # gradient = computePolarization(gradient,
   #                                coefficientField,
   #                                referenceCoefficient
   #                               )
   # gradient = transform(gradient)
   # gradient = applyGammaHat(gradient,
   #                          gamma,
   #                          referenceStiffness,
   #                          lattice
   #                         )
   # gradient = setAverage!(gradient,
   #                        macroscopicGradient,
   #                       )
   # gradient = inverseTransform(gradient)
  end


end

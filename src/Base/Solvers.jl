export solve!

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
          approximationMethod.ansatzSpace,
          problem.lattice
         )
  postprocess!(problem)
  problem
end

export ConvergenceCriterion,
       CauchyConvergenceCriterion,
       NormConvergenceCriterion,
       init!,
       set!,
       computeError

abstract type ConvergenceCriterion end

type CauchyConvergenceCriterion <: ConvergenceCriterion
  oldGradient :: Nullable{SolutionTensorField}

  function CauchyConvergenceCriterion()
    new(Nullable{SolutionTensorField}())
  end
end

function init!(
               c :: CauchyConvergenceCriterion,
               f :: SolutionTensorField{R,F,N}
              ) where {R,F,N}
  c.oldGradient = Nullable(copy(f))
  c
end

function set!(
              c :: CauchyConvergenceCriterion,
              f :: SolutionTensorField{R,F,N}
             ) where {R,F,N}
  copy!(get(c.oldGradient),f)
  c
end

function computeError(
                      c :: CauchyConvergenceCriterion,
                      f :: SolutionTensorField{R,F,N}
                     ) where {R,F,N}
  norm(get(c.oldGradient)-f)/norm(f)
end

type NormConvergenceCriterion <: ConvergenceCriterion
  oldNorm :: Float64

  function NormConvergenceCriterion()
    new(0.0)
  end
end

function init!(
               c :: NormConvergenceCriterion,
               f :: SolutionTensorField{R,F,N}
              ) where {R,F,N}
  c.oldNorm = norm(f)
  c
end

function set!(
              c :: NormConvergenceCriterion,
              f :: SolutionTensorField{R,F,N}
             ) where {R,F,N}
  c.oldNorm = norm(f)
  c
end


function computeError(
                      c :: NormConvergenceCriterion,
                      f :: SolutionTensorField{R,F,N}
                     ) where {R,F,N}
  n = norm(f)
  norm(c.oldNorm-n)/n
end

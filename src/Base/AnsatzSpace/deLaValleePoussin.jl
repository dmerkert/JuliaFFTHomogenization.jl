export deLaValleePoussinMeansSpace

type deLaValleePoussinMeansSpace{R <: AbstractFloat} <: OrthonormalTranslationInvariantSpace
  slope :: Array{R,1}

  function deLaValleePoussinMeansSpace(slope :: Array{R,1}) where {R}
    @argcheck all(slope .>= 0.0)
    @argcheck all(slope .<= 0.5)

    new{R}(slope)
  end
end

function bracketSum(
                    frequency :: Array{I,1},
                    space :: deLaValleePoussinMeansSpace{R},
                    L :: Lattice
                   ) where {
                            I <: Integer,
                            R
                           }
  1.0
end

function ck(
            frequency :: Array{I,1},
            space :: deLaValleePoussinMeansSpace{R},
            L :: Lattice
           ) where {
                    I <: Integer,
                    R
                   }
  delaValleePoussinMean(
                        frequency,
                        L,
                        space.slope,
                        true
                       )
end

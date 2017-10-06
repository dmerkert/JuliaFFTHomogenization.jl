export BoxSplineSpace, ck, bracketSumRange

type BoxSplineSpace{N} <: OrthonormalTranslationInvariantSpace
  Xi :: NTuple{N,Array{Float64,1}}
  bSRange :: Tuple{Int,Int}

  function BoxSplineSpace(Xi :: NTuple{N,Array{Float64,1}},bSRange ::
                          Tuple{Int,Int} = (-5,5)) where {N}
    new{N}(Xi,bSRange)
  end
end

function ck(
            frequency :: Array{I,1},
            space :: BoxSplineSpace{N},
            L :: Lattice
           ) where {
                    I <: Integer,
                    N
                   }


  ckPeriodicBoxSpline(
                      L.MTFactorize\frequency,
                      space.Xi
                     )
end

function bracketSumRange(space :: BoxSplineSpace)
  space.bSRange
end

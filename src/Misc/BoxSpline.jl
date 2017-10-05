export ckPeriodicBoxSpline

function ckPeriodicBoxSpline(
                             k :: Array{Float64,1},
                             Xi :: NTuple{N,Array{Float64,1}}
                            ) where {N}
  prod(
       broadcast(
                 x -> sinc.(x'*k),
                 Xi
                )
      )
end

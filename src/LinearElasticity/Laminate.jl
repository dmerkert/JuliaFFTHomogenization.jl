export laminateFormula

#= function laminateFormula( =#
#=                             C1 :: T1, =#
#=                             C2 :: T2, =#
#=                             amountC1 :: R, =#
#=                             amountC2 :: R, =#
#=                             n :: Array{R,1} =#
#=                            ) where { =#
#=                                     T1 <: StiffnessTensor, =#
#=                                     T2 <: StiffnessTensor, =#
#=                                     R <: AbstractFloat =#
#=                                    } =#

#=   amountC1 /= (amountC1+amountC2) =#
#=   amountC2 /= (amountC1+amountC2) =#

#=   n /= norm(n) =#

#=   P = zeros(6,6) =#

#=   toInverseVoigtSSymFourthOrder!( =#
#=                                  P, =#
#=                                  (i,j,k,l) -> =#
#=                                  _evaluateProjection(n,i,j,k,l) =#
#=                                 ) =#

#=   Id = eye(6) =#
#=   Id[4,4] = 0.5 =#
#=   Id[5,5] = 0.5 =#
#=   Id[6,6] = 0.5 =#

#=   eigC1 = eig(C1) =#
#=   eigC2 = eig(C2) =#

#=   λ = max(max(eigC1...),max(eigC2...)) + 1.0 =#

#=   CEffective = =#
#=   amountC1* =#
#=     inv( =#
#=         P + λ*inv(convert(AnisotropicStiffnessTensor,C1).C - λ*Id) =#
#=        ) + =#
#=   amountC2* =#
#=     inv( =#
#=         P + λ*inv(convert(AnisotropicStiffnessTensor,C2).C - λ*Id) =#
#=     ) =#

#=   CEffective = inv(CEffective) =#
#=   CEffective -= P =#
#=   CEffective ./= λ =#
#=   CEffective = inv(CEffective) =#
#=   CEffective += λ*Id =#
#=   AnisotropicStiffnessTensor(CEffective) =#
#= end =#

#= function _evaluateProjection( =#
#=                               n :: Array{R,1}, =#
#=                               i :: I, =#
#=                               j :: I, =#
#=                               l :: I, =#
#=                               m :: I =#
#=                              ) where { =#
#=                                       R <: AbstractFloat, =#
#=                                       I <: Integer =#
#=                                      } =#
#=   0.5( =#
#=       (j==l)*n[i]*n[m] + =# 
#=       (j==m)*n[i]*n[l] + =# 
#=       (i==l)*n[j]*n[m] + =# 
#=       (i==m)*n[j]*n[l] =#
#=      ) - n[i]*n[j]*n[l]*n[m] =#
#= end =#

function laminateFormula(
                         Stiffnesses :: Array{S,1},
                         volumes :: Array{R,1},
                         normal :: Array{R,1}
                        ) where {S,R}

  P = Array{R}(6,6)
  toInverseVoigtSSymFourthOrder!(
                                 P,
                                 (i,j,l,m) ->
                                 0.5(
                                     (j==l) *
                                     normal[i] *
                                     normal[m] +
                                     (j==m) *
                                     normal[i] *
                                     normal[l] +
                                     (i==l) *
                                     normal[j] *
                                     normal[m] +
                                     (i==m) *
                                     normal[j] *
                                     normal[l]
                                    ) -
                                 normal[i] *
                                 normal[j] *
                                 normal[l] *
                                 normal[m]
                                )
  #= @assert det(P) != 0.0 =#

  PC = AnisotropicStiffnessTensor(P)

  eigenvalues = eig.(Stiffnesses)
  λ = -0.5(min(min.(eigenvalues...)...) + max(max.(eigenvalues...)...))
  #= λ = 2.0max(max.(eigenvalues...)...)+1.0 =#

  λId = λ*IsotropicStiffnessTensor(LamesFirstParameter(0.0),ShearModulus(0.5))

  #= C = Stiffnesses .- λId =#
  #= C = inv.(C) =#
  #= C = PC .+ λ*C =#
  #= C = inv.(C) =#
  #= C = volumes .* C =#
  #= C = sum(C) =#
  #= C = inv(C) =#
  #= C = C - PC =#
  #= C = 1.0/λ * C =#
  #= C = inv(C) =#
  #= C = C + λId =# 

  A = (
       inv(
           1.0/λ .* (
                     inv(
                         sum(
                             volumes.*
                             inv.(
                                  PC .+ λ*inv.(Stiffnesses .- λId)
                                 )
                            )
                        ) .- PC
                    )
          ) .+ λId 
      ) :: AnisotropicStiffnessTensor

  0.5(A+A')

  #= inv( =#
  #=     1.0/λ .* ( =#
  #=               inv( =#
  #=                   sum( =#
  #=                       volumes.* =#
  #=                       inv.( =#
  #=                            PC .+ λ*inv.(Stiffnesses .- λId) =#
  #=                           ) =#
  #=                      ) =#
  #=                  ) .- PC =#
  #=              ) =#
  #=    ) .+ λId :: AnisotropicStiffnessTensor =#
end

function _checkStiff(Stiffnesses)
  yyy = 0.0
  for i in Stiffnesses
    yyy += norm(convert(AnisotropicStiffnessTensor,i).C)
  end

  yyy
end

function symmetricity(C :: Array{R,2}) where {R}
  norm(C - 0.5(C + C'))/norm(C)
end

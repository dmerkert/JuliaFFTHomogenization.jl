export StiffnessTensor,
       IsotropicStiffnessTensor,
       TransversalIsotropicZStiffnessTensor,
       DiagonalStiffnessTensor,
       AnisotropicStiffnessTensor,
       eig,
       mult!,
       convert,
       CompositeArithmeticMeanStiffnessTensor,
       CompositeHarmonicMeanStiffnessTensor,
       CompositeAvgArithmeticHarmonicMeanTensor,
       CompositeLaminateStiffnessTensor,
       inv,
       *


abstract type StiffnessTensor <: CoefficientTensor end

immutable IsotropicStiffnessTensor <: StiffnessTensor
  lambda :: LamesFirstParameter
  mu     :: ShearModulus
end
IsotropicStiffnessTensor(a :: StiffnessParameter,
                         b :: StiffnessParameter) =
IsotropicStiffnessTensor(convert(LamesFirstParameter,a,b),convert(ShearModulus,a,b))

immutable TransversalIsotropicZStiffnessTensor <: StiffnessTensor
  E_p   :: YoungsModulus
  E_t   :: YoungsModulus
  nu_p  :: PoissonsRatio
  nu_pt :: PoissonsRatio
  nu_tp :: PoissonsRatio
  mu_t  :: ShearModulus
end

immutable DiagonalStiffnessTensor <: StiffnessTensor
  v :: Array{Float64,1}
  function DiagonalStiffnessTensor(v :: Array{Float64,1})
    @argcheck length(v) == 6

    new(v)
  end
end

immutable AnisotropicStiffnessTensor <: StiffnessTensor
  C :: Array{Float64,2}
  function AnisotropicStiffnessTensor(C :: Array{Float64,2})
    @argcheck size(C) == (6,6)

    new(C)
  end
end

immutable CompositeLaminateStiffnessTensor{S <: StiffnessTensor} <: StiffnessTensor
  Stiffnesses :: Array{S,1}
  volumes :: Array{Float64,1}
  normal :: Array{Float64,1}
  transformation :: Array{Int,2}

  function CompositeLaminateStiffnessTensor{S}(Stiffnesses :: Array{S,1},
                                               volumes :: Array{Float64,1},
                                               normal :: Array{Float64,1},
                                               transformation :: Array{Int,2}
                                              ) where {S}
    @argcheck length(Stiffnesses) > 0
    @argcheck length(Stiffnesses) == length(volumes)
    @argcheck length(normal) == 3
    @argcheck all(volumes >= 0.0)
    @argcheck sum(volumes) > 0.0
    @argcheck round(Int,det(transformation)) != 0
    @argcheck size(transformation,1) == size(transformation,2)
    @argcheck size(transformation,1) == 3

    volumes /= sum(volumes)
    normal /= norm(normal)

    new(
        Stiffnesses,
        volumes,
        normal,
        transformation
       )
  end
end

for tensor in [
               :CompositeArithmeticMeanStiffnessTensor,
               :CompositeHarmonicMeanStiffnessTensor,
               :CompositeAvgArithmeticHarmonicMeanTensor
              ]
  @eval begin
    immutable ($tensor){S <: StiffnessTensor,R} <: StiffnessTensor
      Stiffnesses :: Array{S,1}
      volumes :: Array{R,1}

      function ($tensor)(Stiffnesses :: Array{S,1},
                         volumes :: Array{R,1}
                        ) where {S,R<:Real}
        @argcheck length(Stiffnesses) > 0
        @argcheck length(Stiffnesses) == length(volumes)
        @argcheck all(volumes .>= 0.0)
        @argcheck sum(volumes) > 0.0

        volumes /= sum(volumes)

        new{S,R}(
                 Stiffnesses,
                 volumes
                )
      end

    end
  end
end

#TODO: optimize für jeden Typ
function eig(A :: AnisotropicStiffnessTensor)
  C = A.C
  C[:,3:6] .*= sqrt(2.0)
  C[3:6,:] .*= sqrt(2.0)
  eigvals(C)
end
eig{S <: StiffnessTensor}(A :: S) = eig(convert(AnisotropicStiffnessTensor,A))

function mult!{R}(stress::Stress{R},
                  stiffness::IsotropicStiffnessTensor,
                  strain::Strain{R}
                 )
  mu = stiffness.mu
  lambda = stiffness.lambda

  stress[1] = (2.0mu+lambda)*strain[1] + lambda*(strain[2]+strain[3])
  stress[2] = (2.0mu+lambda)*strain[2] + lambda*(strain[1]+strain[3])
  stress[3] = (2.0mu+lambda)*strain[3] + lambda*(strain[1]+strain[2])
  stress[4] = mu*strain[4]
  stress[5] = mu*strain[5]
  stress[6] = mu*strain[6]
  stress
end

function mult!{R}(stress::Stress{R},
                  stiffness::DiagonalStiffnessTensor,
                  strain::Strain{R}
                 )
  stress.val = stiffness.v.*strain.val
  stress
end

function mult!{R}(stress::Stress{R},
                  stiffness::TransversalIsotropicZStiffnessTensor,
                  strain::Strain{R}
                 )
  nu_p = stiffness.nu_p
  nu_pt = stiffness.nu_pt
  nu_tp = stiffness.nu_tp
  E_p = stiffness.E_p
  E_t = stiffness.E_t
  mu_t = stiffness.mu_t

  Delta = (1+nu_p)*(1-nu_p-2nu_pt*nu_tp)/(E_p^2*E_t)

  c11 = (1-nu_pt*nu_tp)/(E_p*E_t*Delta)
  c12 = (nu_p+nu_tp*nu_pt)/(E_p*E_t*Delta)
  c13 = (nu_tp+nu_p*nu_tp)/(E_p*E_t*Delta)
  c31 = (nu_pt+nu_p*nu_pt)/(E_p^2*Delta)
  c32 = (nu_pt*(1+nu_p))/(E_p^2*Delta)
  c33 = (1-nu_p^2)/(E_p^2*Delta)
  c44 = mu_t
  c66 = 0.5E_p/(1+nu_p)

  stress.val[1] = c11*strain.val[1]+c12*strain.val[2]+c13*strain.val[3]
  stress.val[2] = c12*strain.val[1]+c11*strain.val[2]+c13*strain.val[3]
  stress.val[3] = c31*strain.val[1]+c32*strain.val[2]+c33*strain.val[3]
  stress.val[4] = c44*strain.val[4]
  stress.val[5] = c44*strain.val[5]
  stress.val[6] = c66*strain.val[6]
  stress
end

function mult!{R}(stress::Stress{R},
                  stiffness::AnisotropicStiffnessTensor,
                  strain::Strain{R}
                 )
  stress.val = stiffness.C*strain.val
  stress
end

function convert(::Type{AnisotropicStiffnessTensor},
                 stiffness::IsotropicStiffnessTensor
  )
  lambda = stiffness.lambda
  mu = stiffness.mu
  C = zeros(6,6)

  C[1,1] = 2.0mu + lambda
  C[2,2] = 2.0mu + lambda
  C[3,3] = 2.0mu + lambda
  C[1,2] = lambda
  C[1,3] = lambda
  C[2,1] = lambda
  C[2,3] = lambda
  C[3,1] = lambda
  C[3,2] = lambda
  C[4,4] = mu
  C[5,5] = mu
  C[6,6] = mu

  AnisotropicStiffnessTensor(C)
end

function convert(::Type{AnisotropicStiffnessTensor},
                 stiffness::DiagonalStiffnessTensor
  )
  C = diagm(stiffness.v)
  AnisotropicStiffnessTensor(C)
end

function convert(::Type{AnisotropicStiffnessTensor},
                 stiffness::TransversalIsotropicZStiffnessTensor
  )
  nu_p = stiffness.nu_p
  nu_pt = stiffness.nu_pt
  nu_tp = stiffness.nu_tp
  E_p = stiffness.E_p
  E_t = stiffness.E_t
  mu_t = stiffness.mu_t

  Delta = (1+nu_p)*(1-nu_p-2nu_pt*nu_tp)/(E_p^2*E_t)

  c11 = (1-nu_pt*nu_tp)/(E_p*E_t*Delta)
  c12 = (nu_p+nu_tp*nu_pt)/(E_p*E_t*Delta)
  c13 = (nu_tp+nu_p*nu_tp)/(E_p*E_t*Delta)
  c31 = (nu_pt+nu_p*nu_pt)/(E_p^2*Delta)
  c32 = (nu_pt*(1+nu_p))/(E_p^2*Delta)
  c33 = (1-nu_p^2)/(E_p^2*Delta)
  c44 = mu_t
  c66 = 0.5E_p/(1+nu_p)

  C = zeros(6,6)
  C[1,1] = (1-nu_pt*nu_tp)/(E_p*E_t*Delta)
  C[1,2] = (nu_p+nu_tp*nu_pt)/(E_p*E_t*Delta)
  C[1,3] = (nu_tp+nu_p*nu_tp)/(E_p*E_t*Delta)
  C[2,1] = C[1,2]
  C[2,2] = C[1,1]
  C[2,3] = C[1,3]
  C[3,1] = (nu_pt+nu_p*nu_pt)/(E_p^2*Delta)
  C[3,2] = (nu_pt*(1+nu_p))/(E_p^2*Delta)
  C[3,3] = (1-nu_p^2)/(E_p^2*Delta)
  C[4,4] = mu_t
  C[5,5] = mu_t
  C[6,6] = 0.5E_p/(1+nu_p)

  AnisotropicStiffnessTensor(C)
end

function convert(::Type{TransversalIsotropicZStiffnessTensor},
                 stiffness::IsotropicStiffnessTensor
                )
  lambda = stiffness.lambda
  mu = stiffness.mu
  E = convert(YoungsModulus,lambda,mu)
  nu = convert(PoissonsRatio,lambda,mu)

  TransversalIsotropicZStiffnessTensor(E,E,nu,nu,nu,mu)
end

convert(::Type{AnisotropicStiffnessTensor},
        stiffness::CompositeArithmeticMeanStiffnessTensor
       ) =
convert(AnisotropicStiffnessTensor,
        sum(stiffness.volumes.*stiffness.Stiffnesses)
       ) :: AnisotropicStiffnessTensor

convert(::Type{AnisotropicStiffnessTensor},
        stiffness::CompositeHarmonicMeanStiffnessTensor
       ) =
convert(AnisotropicStiffnessTensor,
        inv(
            sum(
                stiffness.volumes.*inv.(stiffness.Stiffnesses)
               )
           )
       ) :: AnisotropicStiffnessTensor


convert(::Type{AnisotropicStiffnessTensor},
        stiffness::CompositeAvgArithmeticHarmonicMeanTensor
       ) =
0.5*(
     sum(stiffness.volumes.*stiffness.Stiffnesses)+
     inv(
         sum(
             stiffness.volumes.*inv.(stiffness.Stiffnesses)
            )
        )
    ) :: AnisotropicStiffnessTensor


promote_rule(::Type{TransversalIsotropicZStiffnessTensor},
             ::Type{IsotropicStiffnessTensor}
            ) = TransversalIsotropicZStiffnessTensor

promote_rule{S <: StiffnessTensor}(::Type{AnisotropicStiffnessTensor},
             ::Type{S}
            ) = AnisotropicStiffnessTensor


for tensor in [
               :CompositeLaminateStiffnessTensor,
               :CompositeArithmeticMeanStiffnessTensor,
               :CompositeHarmonicMeanStiffnessTensor,
               :CompositeAvgArithmeticHarmonicMeanTensor
              ]

  @eval promote_rule{S <: StiffnessTensor}(::Type{($tensor)},
                                           ::Type{S}
                                          ) = AnisotropicStiffnessTensor
end

promote_rule{S <: StiffnessTensor}(::Type{DiagonalStiffnessTensor},
             ::Type{S}
            ) = AnisotropicStiffnessTensor

for i in [:+,:-]
  @eval ($i)(A :: IsotropicStiffnessTensor,
             B :: IsotropicStiffnessTensor) =
  IsotropicStiffnessTensor(($i)(A.lambda,B.lambda),
                           ($i)(A.mu,B.mu))

  @eval ($i)(A :: TransversalIsotropicZStiffnessTensor,
             B :: TransversalIsotropicZStiffnessTensor) =
  TransversalIsotropicZStiffnessTensor(($i)(A.E_p,B.E_p),
                                       ($i)(A.E_t,B.E_t),
                                       ($i)(A.nu_p,B.nu_p),
                                       ($i)(A.nu_pt,B.nu_pt),
                                       ($i)(A.nu_tp,B.nu_tp),
                                       ($i)(A.mu_t,B.mu_t) 
                                      )

  @eval ($i)(A :: DiagonalStiffnessTensor,
             B :: DiagonalStiffnessTensor) =
  DiagonalStiffnessTensor(($i)(A.v,B.v))

  @eval ($i)(A :: AnisotropicStiffnessTensor,
             B :: AnisotropicStiffnessTensor) =
  AnisotropicStiffnessTensor(($i)(A.C,B.C))

  @eval ($i)(A :: StiffnessTensor,
             B :: StiffnessTensor) =
  ($i)(promote(A,B)...)
end

*(a :: R, A :: IsotropicStiffnessTensor) where {R <: Real} =
IsotropicStiffnessTensor(
                         LamesFirstParameter(a*A.lambda),
                         ShearModulus(a*A.mu)
                        )

*(a :: R, A :: DiagonalStiffnessTensor) where {R <: Real} =
DiagonalStiffnessTensor(a*A.v)

*(a :: R, A :: S) where {R <: Real,S <: StiffnessTensor} =
AnisotropicStiffnessTensor(a*convert(AnisotropicStiffnessTensor,A).C)


Base.inv(A :: DiagonalStiffnessTensor) = DiagonalStiffnessTensor(1.0./A.v)

Base.inv(A :: S) where {S <: StiffnessTensor} =
AnisotropicStiffnessTensor(inv(convert(AnisotropicStiffnessTensor,A).C))

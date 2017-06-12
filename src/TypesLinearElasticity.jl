abstract StiffnessTensor <: CoefficientTensor
abstract StiffnessParameter


type BulkModulus{R <: AbstractFloat} <: StiffnessParameter
  val :: R
  function BulkModulus{R}(val :: R)
    @argcheck val > 0
    new(val)
  end
end
BulkModulus{R <: AbstractFloat}(val :: R) = BulkModulus{R}(val)
BulkModulus() = BulkModulus(1.0)

type YoungsModulus{R <: AbstractFloat} <: StiffnessParameter
  val :: R
end
YoungsModulus() = YoungsModulus(1.0)

type LamesFirstParameter{R <: AbstractFloat} <: StiffnessParameter
  val :: R
  function LamesFirstParameter{R}(val :: R)
    val >= 0.0 || warning("Lame's first parameter is usually positive")
    new(val)
  end
end
LamesFirstParameter{R <: AbstractFloat}(val :: R) = LamesFirstParameter{R}(val)
LamesFirstParameter() = LamesFirstParameter(1.0)

type ShearModulus{R <: AbstractFloat} <: StiffnessParameter
  val :: R
  function ShearModulus{R}(val :: R)
    @argcheck val >= 0.0
    new(val)
  end
end
ShearModulus{R <: AbstractFloat}(val :: R) = ShearModulus{R}(val)
ShearModulus() = ShearModulus(1.0)

type PoissonsRatio{R <: AbstractFloat} <: StiffnessParameter
  val :: R
  function PoissonsRatio{R}(val :: R)
    @argcheck -1.0 <= val <= 0.5
    new(val)
  end
end
PoissonsRatio{R <: AbstractFloat}(val :: R) = PoissonsRatio{R}(val)
PoissonsRatio() = PoissonsRatio(0.25)

function convert!{R <: StiffnessParameter}(a :: R,b :: R,c :: StiffnessParameter)
  a.val = b.val
end
function convert!{R <: StiffnessParameter}(a :: R,b :: StiffnessParameter,c :: R)
  a.val = c.val
end
function convert!(K :: BulkModulus,E :: YoungsModulus,l :: LamesFirstParameter)
  K.val = (E.val+3.0l.val+sqrt(E.val^2+9.0l.val^2+2.0E.val*l.val))/6.0
end
function convert!(K :: BulkModulus,E :: YoungsModulus,mu :: ShearModulus)
  K.val = E.val*mu.val/(3.0*(3.0mu.val-E.val))
end
function convert!(K :: BulkModulus,E :: YoungsModulus, nu :: PoissonsRatio)
  K.val = E.val/(3.0*(1.0-2.0nu.val))
end
function convert!(mu :: ShearModulus, E :: YoungsModulus, nu :: PoissonsRatio)
  mu.val = E.val/(2.0*(1.0+nu.val))
end
function convert!(l :: LamesFirstParameter, K :: BulkModulus, mu :: ShearModulus)
  l.val = K.val-2.0/3.0*mu.val;
end
function convert!(E :: YoungsModulus, K :: BulkModulus, nu :: PoissonsRatio)
  E.val = 2.0*K.val*(1.0-2.0*nu.val);
end




type Strain{R <: Number} <: GradientSolutionTensor
  val :: Array{R,1}

  function Strain{R}(val :: Array{R,1})
    @argcheck length(val) == 6
    new(val)
  end
end
Strain{R}(val :: Array{R,1}) = Strain{R}(val)

function Strain(R)
  Strain(Array(R,6))
end

function Strain{R <: Number}(val :: Array{R,2})
  @argcheck size(val) == (3,3)
  @argcheck issymmetric(val)

  Strain([val[[1,5,9]];2*val[[8,7,4]]])
end

function Base.getindex(strain::Strain, i::Int)
  strain.val[i]
end

function Base.getindex(strain::Strain, I::Vararg{Int, 2})
  @argcheck 1 <= I[1] <= 3
  @argcheck 1 <= I[2] <= 3
  I[1] == I[2] && return strain.val[I[1]]
  I[1] == 1 && I[2] == 2 && return 0.5strain.val[6]
  I[1] == 1 && I[2] == 3 && return 0.5strain.val[5]
  I[1] == 2 && I[2] == 1 && return 0.5strain.val[6]
  I[1] == 2 && I[2] == 3 && return 0.5strain.val[4]
  I[1] == 3 && I[2] == 1 && return 0.5strain.val[5]
  I[1] == 3 && I[2] == 2 && return 0.5strain.val[4]
end

function Base.setindex!{R}(strain::Strain{R}, v::R, i::Int)
  strain.val[i] = v
  strain
end

function Base.setindex!{R}(strain::Strain{R}, v::R, I::Vararg{Int, 2})
  @argcheck 1 <= I[1] <= 3
  @argcheck 1 <= I[2] <= 3
  I[1] == I[2] && (strain.val[I[1]] = v)
  I[1] == 1 && I[2] == 2 && (strain.val[6] = 2.0v)
  I[1] == 1 && I[2] == 3 && (strain.val[5] = 2.0v)
  I[1] == 2 && I[2] == 1 && (strain.val[6] = 2.0v)
  I[1] == 2 && I[2] == 3 && (strain.val[4] = 2.0v)
  I[1] == 3 && I[2] == 1 && (strain.val[5] = 2.0v)
  I[1] == 3 && I[2] == 2 && (strain.val[4] = 2.0v)
  strain
end

function fromVoigt!{R <: Number}(A :: Array{R,2},strain::Strain{R})
  @argcheck size(A) == (3,3)

  A[1,1] = strain.val[1]
  A[2,2] = strain.val[2]
  A[3,3] = strain.val[3]
  A[2,3] = 0.5strain.val[4]
  A[1,3] = 0.5strain.val[5]
  A[1,2] = 0.5strain.val[6]
  A[2,1] = A[1,2]
  A[3,1] = A[1,3]
  A[3,2] = A[2,3]
  A
end

function toVoigt!{R <: Number}(strain :: Strain{R}, A :: Array{R,2})
  @argcheck size(A) == (3,3)

  strain.val[1] = A[1,1]
  strain.val[2] = A[2,2]
  strain.val[3] = A[3,3]
  strain.val[4] = 2.0A[2,3]
  strain.val[5] = 2.0A[1,3]
  strain.val[6] = 2.0A[1,2]
  strain
end

type Stress{R <: Number} <: FluxSolutionTensor
  val :: Array{R,1}

  function Stress{R}(val :: Array{R,1})
    @argcheck length(val) == 6
    new(val)
  end

end
Stress{R}(val :: Array{R,1}) = Stress{R}(val)

function Stress(R)
  Stress(Array(R,6))
end

function Stress{R <: Number}(val :: Array{R,2})
  @argcheck size(val) == (3,3)
  @argcheck issymmetric(val)

  Stress([val[[1,5,9]];val[[8,7,4]]])
end

function Base.getindex(stress::Stress, i::Int)
  stress.val[i]
end

function Base.getindex(stress::Stress, I::Vararg{Int, 2})
  @argcheck 1 <= I[1] <= 3
  @argcheck 1 <= I[2] <= 3
  I[1] == I[2] && return stress.val[I[1]]
  I[1] == 1 && I[2] == 2 && return stress.val[6]
  I[1] == 1 && I[2] == 3 && return stress.val[5]
  I[1] == 2 && I[2] == 1 && return stress.val[6]
  I[1] == 2 && I[2] == 3 && return stress.val[4]
  I[1] == 3 && I[2] == 1 && return stress.val[5]
  I[1] == 3 && I[2] == 2 && return stress.val[4]
end

function Base.setindex!{R}(stress::Stress{R}, v::R, i::Int)
  stress.val[i] = v
  stress
end

function Base.setindex!{R}(stress::Stress{R}, v::R, I::Vararg{Int, 2})
  @argcheck 1 <= I[1] <= 3
  @argcheck 1 <= I[2] <= 3
  I[1] == I[2] && (stress.val[I[1]] = v)
  I[1] == 1 && I[2] == 2 && (stress.val[6] = v)
  I[1] == 1 && I[2] == 3 && (stress.val[5] = v)
  I[1] == 2 && I[2] == 1 && (stress.val[6] = v)
  I[1] == 2 && I[2] == 3 && (stress.val[4] = v)
  I[1] == 3 && I[2] == 1 && (stress.val[5] = v)
  I[1] == 3 && I[2] == 2 && (stress.val[4] = v)
  stress
end

function fromVoigt!{R <: Number}(A :: Array{R,2},stress::Stress{R})
  @argcheck size(A) == (3,3)

  A[1,1] = stress.val[1]
  A[2,2] = stress.val[2]
  A[3,3] = stress.val[3]
  A[2,3] = stress.val[4]
  A[1,3] = stress.val[5]
  A[1,2] = stress.val[6]
  A[2,1] = A[1,2]
  A[3,1] = A[1,3]
  A[3,2] = A[2,3]
  A
end

function toVoigt!{R <: Number}(stress :: Stress{R}, A :: Array{R,2})
  @argcheck size(A) == (3,3)

  stress.val[1] = A[1,1]
  stress.val[2] = A[2,2]
  stress.val[3] = A[3,3]
  stress.val[4] = A[2,3]
  stress.val[5] = A[1,3]
  stress.val[6] = A[1,2]
  stress
end

type Displacement{R <: Number} <: PrimarySolutionTensor
  val :: Array{R,1}

  function Displacement{R}(val :: Array{R,1})
    @argcheck length(val) == 3
    new(val)
  end
end
Displacement{R}(val :: Array{R,1}) = Displacement{R}(val)

function Displacement(R)
  Displacement(Array(R,3))
end

immutable IsotropicStiffnessTensor{R <: AbstractFloat} <: StiffnessTensor
  lambda :: LamesFirstParameter{R}
  mu     :: ShearModulus{R}
end

immutable TransversalIsotropicZ{R <: AbstractFloat} <: StiffnessTensor
  E_p   :: YoungsModulus{R}
  E_t   :: YoungsModulus{R}
  nu_p  :: PoissonsRatio{R}
  nu_pt :: PoissonsRatio{R}
  nu_tp :: PoissonsRatio{R}
  mu_t  :: ShearModulus{R}
end

immutable DiagonalStiffnessTensor{R <: AbstractFloat} <: StiffnessTensor
  v :: Array{R,1}
end

type StiffnessTensorField{R <: StiffnessTensor} <: CoefficientTensorField
  val :: Array{R}
  function StiffnessTensorField{R}(val:: Array{R})
    new(val)
  end
end
StiffnessTensorField(dims) = 
StiffnessTensorField{StiffnessTensor}(Array{StiffnessTensor}(dims))
StiffnessTensorField(R,dims) = 
StiffnessTensorField{R}(Array{R}(dims))
get(stiffnessTensorField :: StiffnessTensorField,Index) =
stiffnessTensorField.val[Index]
size(stiffnessTensorField :: StiffnessTensorField) =
size(stiffnessTensorField.val)
function get!{S <: StiffnessTensor}(stiffness :: S,stiffnessTensorField ::
                 StiffnessTensorField,Index)
  stiffness = stiffnessTensorField.val[Index,:]
end
function set!{S <: StiffnessTensor}(stiffnessTensorField :: StiffnessTensorField,stiffness :: S,Index)
  stiffnessTensorField.val[Index,:] = stiffness
end

type StrainField{R <: Number} <: GradientSolutionTensorField
  val :: Array{R}
end
StrainField(R,dims) = StrainField(Array{R}(tuple(dims...,6)))
function get!{R}(strain :: Strain{R},strainField :: StrainField{R},Index)
  strain.val = strainField.val[Index,:]
end
function set!{R}(strainField :: StrainField{R},strain :: Strain{R},Index)
  strainField.val[Index,:] = strain.val
end
function add!{R}(strain :: Strain{R},strainField :: StrainField{R},Index)
  strain.val += strainField.val[Index,:]
end
function Base.setindex!{R}(strainField::StrainField{R}, strain::Strain{R},
                           I::CartesianIndex)
  set!(strainField,strain,I)
end
size(strainField :: StrainField) = size(strainField.val)[1:end-1]

type StressField{R <: Number} <: FluxSolutionTensorField
  val :: Array{R}
end
StressField(R,dims) = StressField(Array{R}(tuple(dims...,6)))
function get!{R}(stress :: Stress{R},stressField :: StressField{R},Index)
  stress.val = stressField.val[Index,:]
end
function set!{R}(stressField :: StressField{R},stress :: Stress{R},Index)
  stressField.val[Index,:] = stress.val
end
function add!{R}(stress :: Stress{R},stressField :: StressField{R},Index)
  stress.val += stressField.val[Index,:]
end
size(stressField :: StressField) = size(stressField.val)[1:end-1]


type DisplacementField{R <: Number} <: PrimarySolutionTensorField
  val :: Array{R}
end
DisplacementField(R,dims) =
DisplacementField(Array{R}(tuple(dims...,6)))


function Transform!{R}(solutionTensor :: SolutionTensor, A :: Array{R,2})
  @argcheck size(A) == (3,3)
  @argcheck typeof(solutionTensor.val[1]) == R

  solutionTensorMatrix = zeros(A)

  fromVoigt!(solutionTensorMatrix,solutionTensor)
  solutionTensorMatrix = A*solutionTensorMatrix*A'
  toVoigt!(solutionTensor,solutionTensorMatrix)
end

__precompile__()
module JuliaFFTHomogenization

import Base.*
import Base./
import Base.+
import Base.-
import Base.==
import Base.getindex
import Base.norm
import Base.promote_rule
import Base.convert
import Base.length
import Base.copy
import Base.copy!
import Base.setindex!
import Base.getindex
using ArgCheck
using MPAWL
import Base.size
using Polynomials

FFTW.set_num_threads(Sys.CPU_CORES)

include("Base/Types.jl")
include("Base/AnsatzSpace.jl")
include("Base/AnsatzSpace/TruncatedTrigonometricPolynomials.jl")
include("Base/AnsatzSpace/translationInvariantSpace.jl")
include("Base/AnsatzSpace/deLaValleePoussin.jl")
include("Base/AnsatzSpace/boxSplines.jl")
include("Base/ApproximationMethod.jl")
include("Base/ConvergenceCriterion.jl")
include("Base/Solvers.jl")
include("Base/Solvers/BasicScheme.jl")

include("LinearElasticity/StiffnessParameters.jl")
include("LinearElasticity/StrainStressDisplacement.jl")
include("LinearElasticity/StrainStressDisplacementFields.jl")
include("LinearElasticity/StiffnessTensors.jl")
include("LinearElasticity/GreenOperators/Gamma0Elasticity.jl")
include("LinearElasticity/Problem.jl")
include("LinearElasticity/ReferenceTensor.jl")
include("LinearElasticity/Tensors.jl")

include("LinearElasticity/Laminate.jl")
include("LinearElasticity/Subsampling.jl")

include("Misc/Misc.jl")
include("Misc/Ellipsoids.jl")
include("Misc/BoxSpline.jl")
include("Misc/PeriodicSmoothSpline.jl")

include("ExampleGeometries/ElasticityHashinEllipsoid.jl")
include("ExampleGeometries/Elasticity1DLaminate.jl")
include("ExampleGeometries/ElasticitySingleFiber.jl")
include("ExampleGeometries/ElasticityBoxSpline.jl")
include("ExampleGeometries/ElasticitySpline.jl")
include("ExampleGeometries/ElasticityPeriodicSmoothSpline.jl")

end # module

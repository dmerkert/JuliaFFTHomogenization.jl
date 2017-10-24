using Polynomials

export getPolys,evalPolys


function _coeffPolys(deg,diff,x)
  if diff < 0
    return [0 for n=0:deg]
  else
    return [diff > n ? zero(x) : factorial(diff)*binomial(n,diff)*x^(n-diff) for n=0:deg]
  end
end

function _testPolys(o,p0,p1)
    l = -0.25
    r = 0.25

    t0 = p0
    t1 = p1

    for i=0:o
      @assert isapprox(polyval(t0,l),0.0,atol=1e-14)
      @assert isapprox(polyval(t1,r),0.0,atol=1e-14)
      @assert isapprox(polyval(t0,0.0),polyval(t1,0.0),atol=1e-14)

      t0 = polyder(t0)
      t1 = polyder(t1)
    end

    @assert abs(polyval(t0,l)) > 1.0
    @assert abs(polyval(t1,l)) > 1.0
    @assert abs(polyval(t0,0.0) -  polyval(t1,0.0)) > 1.0
end

function getPolys(o)
  if o == -1
    return (Poly([1]),Poly([1]))
  elseif o == 0
    l = -0.25
    r = 0.25
    m = 100
    A = [
         _coeffPolys(1,0,l)'  _coeffPolys(1,m,r)';
         _coeffPolys(1,m,l)'  _coeffPolys(1,0,r)';
         _coeffPolys(1,0,0)' -_coeffPolys(1,0,0)';
         _coeffPolys(1,0,0)'  _coeffPolys(1,m,0)';
        ]

    b=[zeros(3);1]
    #= return (Poly([0.25;1]),Poly([0.25;-1])) =#

    #= c = round.(Int,A\b) =#
    c = Matrix{Rational{BigInt}}(A)\Vector{Rational{BigInt}}(b)

    p0 = Poly(c[1:2])
    p1 = Poly(c[3:4])

    _testPolys(o,p0,p1)

    return (p0,p1)
  elseif o == 1
    l = -0.25
    r = 0.25
    m = 100
    A = [
         _coeffPolys(3,0,l)'  _coeffPolys(3,m,r)';
         _coeffPolys(3,1,l)'  _coeffPolys(3,m,r)';
         _coeffPolys(3,m,l)'  _coeffPolys(3,0,r)';
         _coeffPolys(3,m,l)'  _coeffPolys(3,1,r)';
         _coeffPolys(3,0,0)' -_coeffPolys(3,0,0)';
         _coeffPolys(3,1,0)' -_coeffPolys(3,1,0)';
         _coeffPolys(3,0,0)'  _coeffPolys(3,m,0)';
         _coeffPolys(3,1,0)'  _coeffPolys(3,m,0)'
        ]
    b=[zeros(6);1;1]

    #= c = round.(Int,A\b) =#
    c = Matrix{Rational{BigInt}}(A)\Vector{Rational{BigInt}}(b)

    p0 = Poly(c[1:4])
    p1 = Poly(c[5:8])

    _testPolys(o,p0,p1)
    return (p0,p1)
  elseif o == 2
    l = -0.25
    r = 0.25
    m = 100
    A = [
         _coeffPolys(4,0,l)'  _coeffPolys(4,m,r)';
         _coeffPolys(4,1,l)'  _coeffPolys(4,m,r)';
         _coeffPolys(4,2,l)'  _coeffPolys(4,m,r)';
         _coeffPolys(4,m,l)'  _coeffPolys(4,0,r)';
         _coeffPolys(4,m,l)'  _coeffPolys(4,1,r)';
         _coeffPolys(4,m,l)'  _coeffPolys(4,2,r)';
         _coeffPolys(4,0,0)' -_coeffPolys(4,0,0)';
         _coeffPolys(4,1,0)' -_coeffPolys(4,1,0)';
         _coeffPolys(4,2,0)' -_coeffPolys(4,2,0)';
         _coeffPolys(4,0,0)'  _coeffPolys(4,m,0)'
        ]

    b=[zeros(9);1]

    #= c = round.(Int,A\b) =#
    c = Matrix{Rational{BigInt}}(A)\Vector{Rational{BigInt}}(b)


    p0 = Poly(c[1:5])
    p1 = Poly(c[6:10])

    _testPolys(o,p0,p1)
    return (p0,p1)
  elseif o == 3
    l = -0.25
    r = 0.25
    m = 100
    A = [
         _coeffPolys(6,0,l)'  _coeffPolys(6,m,r)';
         _coeffPolys(6,1,l)'  _coeffPolys(6,m,r)';
         _coeffPolys(6,2,l)'  _coeffPolys(6,m,r)';
         _coeffPolys(6,3,l)'  _coeffPolys(6,m,r)';
         _coeffPolys(6,m,l)'  _coeffPolys(6,0,r)';
         _coeffPolys(6,m,l)'  _coeffPolys(6,1,r)';
         _coeffPolys(6,m,l)'  _coeffPolys(6,2,r)';
         _coeffPolys(6,m,l)'  _coeffPolys(6,3,r)';
         _coeffPolys(6,0,0)' -_coeffPolys(6,0,0)';
         _coeffPolys(6,1,0)' -_coeffPolys(6,1,0)';
         _coeffPolys(6,2,0)' -_coeffPolys(6,2,0)';
         _coeffPolys(6,3,0)' -_coeffPolys(6,3,0)';
         _coeffPolys(6,0,0)'  _coeffPolys(6,m,0)';
         _coeffPolys(6,1,0)'  _coeffPolys(6,m,0)'
        ]

    b = [zeros(12);1;2]

    #= c = round.(Int,A\b) =#
    c = Matrix{Rational{BigInt}}(A)\Vector{Rational{BigInt}}(b)


    p0 = Poly(c[1:7])
    p1 = Poly(c[8:14])

    _testPolys(o,p0,p1)
    return (p0,p1)
  elseif o == 4
    l = -0.25
    r = 0.25
    m = 100
    A = [
         _coeffPolys(7,0,l)'  _coeffPolys(7,m,r)';
         _coeffPolys(7,1,l)'  _coeffPolys(7,m,r)';
         _coeffPolys(7,2,l)'  _coeffPolys(7,m,r)';
         _coeffPolys(7,3,l)'  _coeffPolys(7,m,r)';
         _coeffPolys(7,4,l)'  _coeffPolys(7,m,r)';
         _coeffPolys(7,m,l)'  _coeffPolys(7,0,r)';
         _coeffPolys(7,m,l)'  _coeffPolys(7,1,r)';
         _coeffPolys(7,m,l)'  _coeffPolys(7,2,r)';
         _coeffPolys(7,m,l)'  _coeffPolys(7,3,r)';
         _coeffPolys(7,m,l)'  _coeffPolys(7,4,r)';
         _coeffPolys(7,0,0)' -_coeffPolys(7,0,0)';
         _coeffPolys(7,1,0)' -_coeffPolys(7,1,0)';
         _coeffPolys(7,2,0)' -_coeffPolys(7,2,0)';
         _coeffPolys(7,3,0)' -_coeffPolys(7,3,0)';
         _coeffPolys(7,4,0)' -_coeffPolys(7,4,0)';
         _coeffPolys(7,0,0)'  _coeffPolys(7,m,0)';
        ]

    b = [zeros(15);1]

    c = Matrix{Rational{BigInt}}(A)\Vector{Rational{BigInt}}(b)

    p0 = Poly(c[1:8])
    p1 = Poly(c[9:16])

    _testPolys(o,p0,p1)
    return (p0,p1)
  elseif o == 5
    l = -0.25
    r = 0.25
    m = 100
    A = [
         _coeffPolys(9,0,l)'  _coeffPolys(9,m,r)';
         _coeffPolys(9,1,l)'  _coeffPolys(9,m,r)';
         _coeffPolys(9,2,l)'  _coeffPolys(9,m,r)';
         _coeffPolys(9,3,l)'  _coeffPolys(9,m,r)';
         _coeffPolys(9,4,l)'  _coeffPolys(9,m,r)';
         _coeffPolys(9,5,l)'  _coeffPolys(9,m,r)';
         _coeffPolys(9,m,l)'  _coeffPolys(9,0,r)';
         _coeffPolys(9,m,l)'  _coeffPolys(9,1,r)';
         _coeffPolys(9,m,l)'  _coeffPolys(9,2,r)';
         _coeffPolys(9,m,l)'  _coeffPolys(9,3,r)';
         _coeffPolys(9,m,l)'  _coeffPolys(9,4,r)';
         _coeffPolys(9,m,l)'  _coeffPolys(9,5,r)';
         _coeffPolys(9,0,0)' -_coeffPolys(9,0,0)';
         _coeffPolys(9,1,0)' -_coeffPolys(9,1,0)';
         _coeffPolys(9,2,0)' -_coeffPolys(9,2,0)';
         _coeffPolys(9,3,0)' -_coeffPolys(9,3,0)';
         _coeffPolys(9,4,0)' -_coeffPolys(9,4,0)';
         _coeffPolys(9,5,0)' -_coeffPolys(9,5,0)';
         _coeffPolys(9,0,0)'  _coeffPolys(9,m,0)';
         _coeffPolys(9,1,0)'  _coeffPolys(9,m,0)';
        ]

    b = [zeros(18);1;1]

    #= c = A\b =#
    c = Matrix{Rational{BigInt}}(A)\Vector{Rational{BigInt}}(b)

    p0 = Poly(c[1:10])
    p1 = Poly(c[11:20])

    _testPolys(o,p0,p1)
    return (p0,p1)
  end
end

function evalPolys(x::Float64,o :: Int)
  (f0,f1) = getPolys(o)
  1.0 + ((x >= -0.25) & (x < 0.0)) * Float64(polyval(f0,x)) + ((x >= 0.0) & (x < 0.25)) * Float64(polyval(f1,x))
end

function evalPolys(x::Float64,f0,f1)
  1.0 + ((x >= -0.25) & (x < 0.0)) * Float64(polyval(f0,x)) + ((x >= 0.0) & (x < 0.25)) * Float64(polyval(f1,x))
end

function evalPolys(x::Vector{Float64},f0,f1)
  prod(evalPolys.(x,f0,f1))
end

function evalPolys(x::Vector{Float64},o :: Int)
  prod(evalPolys.(x,o))
end



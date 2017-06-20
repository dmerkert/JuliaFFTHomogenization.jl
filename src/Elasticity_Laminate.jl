function Elasticity_Laminate{R <: Real,I <: Integer}
  (
   patternMatrix :: Array{I,2};
   Normal :: Array{R,1},
   StripWidth :: R,
   nuStrip :: PoissonsRatio,
   nuGap :: PoissonsRatio,
   EStrip :: YoungsModulus,
   EGap :: YoungsModulus
  )

  
end

function _Elasticity_Laminate!
  (
   Normal :: Array{R,1},
   StripWidth :: R,
   nuStrip :: PoissonsRatio,
   nuGap :: PoissonsRatio,
   EStrip :: YoungsModulus,
   EGap :: YoungsModulus,
   point :: Array{R,1}
  )

  X = point[1]
  Y = point[2]

  normal = Normal/norm(Normal)
  direction = [-normal[2];normal[1]]

  offsetList = PeriodicLineOffsets(normal,'StripPrecision',1e-2*pp.StripWidth)

  %detmine the distance of all points to the lines
  IStrip = false(size(points(1,:)))
  for offsetIndex = 1:size(offsetList,1)
    offset = offsetList(offsetIndex,:)'
    rhs = points(1:2,:)
    rhs(1,:) = rhs(1,:) - offset(1)
    rhs(2,:) = rhs(2,:) - offset(2)

    %solve lambda*direction + b + mu * normal = X
    %i.e. [direction, normal] * [lambdamu] = X - b
    y = [direction,normal]\rhs
    mu = y(2,:)

    %assume norm(normal) == 1
    IStrip = IStrip | ((mu <= 0.5*pp.StripWidth) & (mu > -0.5*pp.StripWidth))
  end

  [lambdaStrip,muStrip] = ENuToLambdaMu(pp.EStrip,pp.nuStrip)
  [lambdaGap,muGap] = ENuToLambdaMu(pp.EGap,pp.nuGap)

  lambda = lambdaGap*ones([N,1])
  lambda(IStrip) = lambdaStrip

  mu = muGap*ones([N,1])
  mu(IStrip) = muStrip

  C = isotropicCLame(lambda,mu,N)
  C = reshape(C,[6,6,1,N])



end

function _PeriodicLineOffsets{R <: Real, I <: Integer}
  (
   normal :: Array{R,1};
   StripPrecision::R=1e-4,
   MaxIterations::I=100,
   DigitsForRounding::I=8
  )

  direction = [-normal[2],normal[1]]

  lineParam = @(lambda,direct,offset) ...
    modM(lambda*direct+offset,eye(2),'Target','symmetric');

offsetList = zeros([1,2]);  %[0,0] is of course a valid offset...

if (direction(1) ~= 0)
    centerLine = @(lambda) lineParam(lambda,direction/direction(1),[0;0]);
    lambdaX = 1;
    while 1==1
        offset = centerLine(lambdaX)';
        for shiftX=-1:1
            for shiftY=-1:1
                offsetList = [offsetList;offset+[shiftX,shiftY];-(offset+[shiftX,shiftY])];
            end
        end
        if norm(centerLine(lambdaX)) < pp.StripPrecision   %full period accieved
            break;
        end
        if lambdaX > pp.MaxIterations
            break;
        end
        lambdaX = lambdaX+1;
    end
end

if (direction(2) ~= 0)
    centerLine = @(lambda) lineParam(lambda,direction/direction(2),[0;0]);
    lambdaY = 1;
    while 1==1
        offset = centerLine(lambdaY)';
        for shiftX=-1:1
            for shiftY=-1:1
                offsetList = [offsetList;offset+[shiftX,shiftY];-(offset+[shiftX,shiftY])];
            end
        end
        
        if norm(centerLine(lambdaY)) < pp.StripPrecision   %full period accieved
            break;
        end
        if lambdaY > pp.MaxIterations
            break;
        end
        lambdaY = lambdaY+1;
    end
end

%account for rounding errors
offsetList = round(offsetList,pp.DigitsForRounding);
offsetList = unique(offsetList,'rows');

end


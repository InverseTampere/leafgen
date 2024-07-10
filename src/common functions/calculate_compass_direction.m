function CD = calculate_compass_direction(CylinderParameters)

% Cylinder start points
startPoints = CylinderParameters.start_point;
% Cylinder axis
cylinderAxis = CylinderParameters.axis;
cylinderAxis = cylinderAxis./sqrt(sum(cylinderAxis.^2,2)); % normalization
% Cylinder length
cylinderLength = CylinderParameters.length;

% Cylinder midpoints
midPoints = startPoints + 0.5*cylinderLength.*cylinderAxis;

% Horizontal projection of midpoints
midPointsHorz = [midPoints(:,1:2) zeros(length(midPoints(:,3)),1)];

% Horizontal mean of cylinder midpoints
meanPoint = mean(midPointsHorz,1);

% Position vectors of horizontal cylinder midpoints wrt mean of midpoints
posVecHorz = midPointsHorz - meanPoint;

% If position vector equals to mean of midpoints, set direction to north
zeroPosInd = (sum(posVecHorz,2)==0);
posVecHorz(:,2) = posVecHorz(:,2) + zeroPosInd;

% Normalized horizontal direction vectors
dirVec = posVecHorz./sqrt(sum(posVecHorz.^2,2));

% Compass directions of cylinder midpoints (y-direction corresponds to
% north)
CD = acos(dirVec*[0 1 0]');
CD(dirVec(:,1)>0) = 2*pi - CD(dirVec(:,1)>0);

end
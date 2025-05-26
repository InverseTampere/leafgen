% This file is part of LeafGen
% 
% LeafGen is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% LeafGen is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with LeafGen.  If not, see <https://www.gnu.org/licenses/>.

function CD = calculate_compass_direction(CylinderParameters)

% Cylinder start points
startPoints = CylinderParameters.start_point;
% Cylinder axis
cylinderAxis = CylinderParameters.axis;
cylinderAxis = cylinderAxis./sqrt(sum(cylinderAxis.^2,2)); % normalization
% Cylinder length
cylinderLength = CylinderParameters.length;
% Maximum cylinder length
maxLen = max(cylinderLength);
% Normalized cylinder lengths
cylLenNorm = cylinderLength./maxLen;

% Cylinder midpoints
midPoints = startPoints + 0.5*cylinderLength.*cylinderAxis;

% Horizontal projection of midpoints
midPointsHorz = [midPoints(:,1:2) zeros(length(midPoints(:,3)),1)];

% Horizontal mean of cylinder midpoints weighted by normalized cylinder 
% length
meanPoint = sum(cylLenNorm.*midPointsHorz,1)./sum(cylLenNorm);

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
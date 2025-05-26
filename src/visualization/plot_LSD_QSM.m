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

function f = plot_LSD_QSM(QSM,Leaves,varargin)

% Initialize values
nBins = 10;
svIntervals = [0 1 0 1 0 2*pi];

% Check additional parameters
i = 1;
NArg = numel(varargin);
while i <= NArg
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'nbins'
                nBins = varargin{i+1};
            case 'variableintervals'
                svIntervals = varargin{i+1};
        end
    end
    i = i + 1;
end

% Initialize figure object
f = figure; clf, hold on

%% Calculate the exact structural variable values for the leaf start points

% Extracting QSM and leaf information
cylinderStartPoint  = QSM.cylinder_start_point;
cylinderAxis        = QSM.cylinder_axis;
cylinderBranchIndex = QSM.cylinder_branch_index;
cylinderLength      = QSM.cylinder_length;
cylinderCount       = QSM.block_count;
petioleStartPoint = Leaves.petiole_start_point;
leafParents       = Leaves.leaf_parent;
leafCount         = Leaves.leaf_count;


% RELATIVE HEIGHT

leafHeights = Leaves.leaf_start_point(:,3);
maxHeight = max([max(leafHeights) max(QSM.cylinder_end_point(:,3)) ...
                 max(QSM.cylinder_start_point(:,3))]);
leafRelHei = leafHeights./maxHeight;


% RELATIVE DISTANCE ALONG SUBBRANCH

leafRelPos = nan(leafCount,1);

% Cylinder index vector
cylIndices = (1:1:cylinderCount)';

% Leaf index vector
leafIndices = (1:1:leafCount)';

for iBranch = 0:max(cylinderBranchIndex)
    % Indices of QSM cylinders belonging in branch
    bcIndices = cylIndices(cylinderBranchIndex == iBranch);
    if sum(bcIndices) == 0
        % the branch has no cylinders
        continue
    end
    % Lengths of corresponding cylinders
    bcLengths = zeros(length(bcIndices),1);
    for j = 1:length(bcLengths)
        bcLengths(j) = cylinderLength(bcIndices(j));
    end
    % Cumulative sum of the lengths
    bcLengthsCumulative = cumsum(bcLengths);

    % Branch cylinder edges relative to subbranch distance
    bcRelEdges = [0; bcLengthsCumulative]/bcLengthsCumulative(end);

    % Loop over branch cylinders
    for iBC = 1:length(bcIndices)
        % Find the indices of the leaves attached to the cylinder
        childLeafInds = leafIndices(leafParents == bcIndices(iBC));
        % If cylinder has no leaves, advance to the next cylinder
        if isempty(childLeafInds) == true
            continue
        end
        % Vetor from cylinder start to petiole start
        cylStartToPetioleStart = petioleStartPoint(childLeafInds,:) ...
                              - cylinderStartPoint(bcIndices(iBC),:);
        % Projection to cylinder axis
        cylAxis = cylinderAxis(bcIndices(iBC),:);
        petioleStartProj = cylStartToPetioleStart*cylAxis' ...
                       .*cylAxis./(norm(cylAxis).^2);
        % Leaf start point values as relative branch distance
        leafRelPos(childLeafInds) = sqrt(sum(petioleStartProj.^2,2)) ...
                                    ./cylinderLength(bcIndices(iBC)) ...
                                    *(bcRelEdges(iBC+1)-bcRelEdges(iBC))...
                                    + bcRelEdges(iBC);
    end
end

% Limit the relative position between 0 and 1
leafRelPos(leafRelPos<0) = 0;
leafRelPos(leafRelPos>1) = 1;

% COMPASS DIRECTION

% Maximum cylinder length
maxLen = max(cylinderLength);
% Normalized cylinder lengths
cylLenNorm = cylinderLength./maxLen;
% Cylinder midpoints
midPoints = cylinderStartPoint + 0.5*cylinderLength.*cylinderAxis;
% Set the mean value of cylinder locations on xy-plane weighted with 
% cylinder length as the origin
xyOrigin = sum(cylLenNorm.*midPoints(:,1:2),1)./sum(cylLenNorm);

% Coordinates of the leaf start points on xy-plane
leafCoordXY = Leaves.leaf_start_point(:,1:2);
% Unit vectors pointing the direction of cylinder start point coordinates
leafDirection = zeros(leafCount,2);
for j = 1:leafCount
    temp = leafCoordXY(j,:) - xyOrigin;
    leafDirection(j,:) = temp/norm(temp);
end
% North set as the zero angle
northDirection = [0 1];
% Angles of the cylinder direction unit vectors w.r.t north
leafComDir = zeros(leafCount,1);
for j = 1:leafCount
    if leafDirection(j,1) <= 0
        leafComDir(j) = acos(dot(northDirection,leafDirection(j,:)));
    else
        leafComDir(j) = 2*pi - acos(dot(northDirection,leafDirection(j,:)));
    end
end

%% Find the sizes of leaves within the given interval

% Indices of leaves on the interval
if all(svIntervals == [0 1 0 1 0 2*pi])
    intervalInds = (1:leafCount)';
else
    h1 = svIntervals(1);
    h2 = svIntervals(2);
    d1 = svIntervals(3);
    d2 = svIntervals(4);
    c1 = svIntervals(5);
    c2 = svIntervals(6);
    boolH = (leafRelHei>=h1) & (leafRelHei<=h2);
    boolD = (leafRelPos>=d1) & (leafRelPos<=d2);
    boolC = (leafComDir>=c1) & (leafComDir<=c2);
    intervalInds = leafIndices(boolH & boolD & boolC);
end

% Calculate the surface areas of individual leaves
leafAreas = Leaves.base_area*Leaves.leaf_scale(intervalInds,1).^2;

% Plot histogram
histogram(leafAreas,linspace(min(leafAreas),max(leafAreas),nBins+1), ...
          'Normalization','pdf','FaceColor',"#00CD94", ...
          'FaceAlpha',0.3)
xlabel("leaf area [m^2]")
ylabel("frequency density")
title("LSD")
axis tight
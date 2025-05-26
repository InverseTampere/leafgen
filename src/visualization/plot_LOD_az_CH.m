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

function f = plot_LOD_az_CH(aShape,Leaves,varargin)

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
            case 'stemcoordinates'
                stemCoordinates = varargin{i+1};
                flagStemCoordinates = true;
        end
    end
    i = i + 1;
end

% Initialize figure object
f = figure; clf, hold on

%% Calculate the exact structural variable values for the leaf start points

% Extracting leaf information
leafCount = Leaves.leaf_count;
leafStartPoints = Leaves.leaf_start_point;
leafNormals = Leaves.leaf_normal;

% RELATIVE HEIGHT

% Heights of the leaves
hLeaves = leafStartPoints(:,3);

% Set the maximum and minimum height based on the point cloud
pCloud = aShape.Points;
h_min = min(pCloud(:,3)); 
h_max = max(pCloud(:,3)); 
hLeaves(hLeaves<h_min) = h_min;
hLeaves(hLeaves>h_max) = h_max;

% Normalize leaf heights
leafRelHei = (hLeaves-h_min)/(h_max-h_min);

% RELATIVE DISTANCE FROM STEM

% Maximum horizontal distance of point cloud from origin
maxHorzDist = max(sqrt(sum(pCloud(:,1:2).^2,2)));

% Initialize variable for relative distance from stem
leafRelDis = zeros(leafCount,1);

% Variables for missed leaves due to problems in defining distance to the
% alphashape edge
nMissed = 0;
indMissed = [];
areaMissed = 0;

if flagStemCoordinates == true % Stem coordinates supplied as input
    for iLeaf = 1:leafCount
        % Index of stem coordinate with z-value below the leaf height
        if stemCoordinates(end,3) < leafStartPoints(iLeaf,3)
            iSC = size(stemCoordinates,1);
        else
            iSC = find(stemCoordinates(:,3) > leafStartPoints(iLeaf,3),1);
        end
        % Relative z-position of leaf between the stem coordinate nodes
        relPos = (leafStartPoints(iLeaf,3)- stemCoordinates(iSC-1,3)) ...
                 /(stemCoordinates(iSC,3) - stemCoordinates(iSC-1,3));
        % Location of the stem center on the height of the leaf
        stemCen = relPos*(stemCoordinates(iSC,:) ...
                          -stemCoordinates(iSC-1,:)) ...
                  + stemCoordinates(iSC-1,:);
        % Horizontal distance from stem to leaf
        stemToLeaf = [leafStartPoints(iLeaf,1:2) 0] - [stemCen(1:2) 0];
        distStemToLeaf = norm(stemToLeaf); 
        % Counter-clockwise angle between north and leaf direction wrt. the
        % stem
        unitStemToLeaf = stemToLeaf/distStemToLeaf;
        if unitStemToLeaf(1) < 0
            compassDir = acos(dot(unitStemToLeaf,[0 1 0]));
        else 
            compassDir = 2*pi - acos(dot(unitStemToLeaf,[0 1 0]));
        end
        % Probing the edge of alpha shape
        nPP = 2*100;
        yCoord = 2*maxHorzDist*linspace(0,1,nPP)';
        initPP = [zeros(nPP,1) yCoord zeros(nPP,1)];
        probePoints = (rotation_matrix([0 0 1],compassDir)*initPP')' ...
                      + stemCen;
        tf = inShape(aShape,probePoints);
        edgeIndex = find(~tf,1,'first') - 1;
        if edgeIndex < 1
            % The stem is outside alphashape
            nMissed = nMissed + 1;
            indMissed = [indMissed, iLeaf];
            areaMissed = areaMissed ...
                         + Leaves.base_area*Leaves.leaf_scale(iLeaf,2).^2;
            continue
        end
        edgeValue = yCoord(edgeIndex);
        % The relative horizontal distance from the stem
        leafRelDis(iLeaf) = distStemToLeaf/edgeValue;
    end
else % Assuming stem to be the z-axis
    for iLeaf = 1:leafCount
        % Unit direction vector on xy-plane
        xyDir = leafStartPoints(iLeaf,1:2) ...
                /sqrt(sum(leafStartPoints(iLeaf,1:2).^2));
        % Probing the edge of alpha shpe
        nPP = 100;
        probeDist = maxHorzDist*linspace(0,1,nPP)';
        probePoints = [probeDist*xyDir, ...
                       leafStartPoints(iLeaf,3)*ones(nPP,1)];
        tf = inShape(aShape,probePoints);
        edgeIndex = find(~tf,1,'first') - 1;
        if edgeIndex < 1
            % The stem is outside alphashape
            nMissed = nMissed + 1;
            indMissed = [indMissed, iLeaf];
            areaMissed = areaMissed ...
                         + Leaves.base_area*Leaves.leaf_scale(iLeaf,2).^2;
            continue
        end
        edgeValue = norm(probePoints(edgeIndex,:));
        % The relative horizontal distance from the stem
        stemToLeaf = [leafStartPoints(iLeaf,1:2) 0];
        leafRelDis(iLeaf) = norm(stemToLeaf)/edgeValue;
    end
end

if nMissed > 0
    warning('%d leaves not included in the \"relative distance from stem\" histogram due to the defined stem being outside alphashape. (Total missed area %.2f m^2)',nMissed,areaMissed)
    leafCount = leafCount - nMissed;
    leafRelDis(indMissed) = nan;
end

% COMPASS DIRECTION

% initialize variable for compass direction wrt. the stem
leafComDir = zeros(leafCount,1);

if flagStemCoordinates == true
    for iLeaf = 1:leafCount
        % Index of stem coordinate with z-value below the leaf height
        if stemCoordinates(end,3) < leafStartPoints(iLeaf,3)
            iSC = size(stemCoordinates,1);
        else
            iSC = find(stemCoordinates(:,3) > leafStartPoints(iLeaf,3),1);
        end
        % Relative z-position of leaf between the stem coordinate nodes
        relPos = (leafStartPoints(iLeaf,3)- stemCoordinates(iSC-1,3)) ...
                 /(stemCoordinates(iSC,3) - stemCoordinates(iSC-1,3));
        % Location of the stem center on the height of the leaf
        stemCen = relPos*(stemCoordinates(iSC,:) ...
                          -stemCoordinates(iSC-1,:)) ...
                  + stemCoordinates(iSC-1,:);
        % Horizontal distance from stem to leaf
        stemToLeaf = [leafStartPoints(iLeaf,1:2) 0] - [stemCen(1:2) 0];
        distStemToLeaf = norm(stemToLeaf); 
        % Counter-clockwise angle between north and leaf direction wrt. the
        % stem
        unitStemToLeaf = stemToLeaf/distStemToLeaf;
        if unitStemToLeaf(1) < 0
            leafComDir(iLeaf) = acos(dot(unitStemToLeaf,[0 1 0]));
        else 
            leafComDir(iLeaf) = 2*pi - acos(dot(unitStemToLeaf,[0 1 0]));
        end
    end
else
    for iLeaf = 1:leafCount
        % Unit horizontal direction vector from the stem to leaf
        stemToLeaf = [leafStartPoints(iLeaf,1:2) 0];
        unitStemToLeaf = stemToLeaf/norm(stemToLeaf);
        % Counter-clockwise compass direction with respect to north 
        % direction
        if unitStemToLeaf(1) < 0
            leafComDir(iLeaf) = acos(dot(unitStemToLeaf,[0 1 0]));
        else 
            leafComDir(iLeaf) = 2*pi - acos(dot(unitStemToLeaf,[0 1 0]));
        end
    end
end

%% Find the azimuth angles of leaves within the given interval

% Leaf index vector
leafIndices = (1:1:leafCount)';

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
    boolD = (leafRelDis>=d1) & (leafRelDis<=d2);
    boolC = (leafComDir>=c1) & (leafComDir<=c2);
    intervalInds = leafIndices(boolH & boolD & boolC);
end

% Calculate azimuth angles of leaf normals within the interval
horzDir = leafNormals(intervalInds,1:2);
normHD  = sqrt(sum(horzDir.^2,2));
horzDir = horzDir((normHD > 1e-6),:);
horzDir = horzDir./sqrt(sum(horzDir.^2,2));
nIntLeaves = size(horzDir,1);
azAngles = nan(nIntLeaves,1);
for iLeaf = 1:nIntLeaves
    if horzDir(iLeaf,1) <= 0
        azAngles(iLeaf) = acos(dot(horzDir(iLeaf,:),[0 1]));
    else
        azAngles(iLeaf) = 2*pi - acos(dot(horzDir(iLeaf,:),[0 1]));
    end
end

% Plot histogram
histogram(azAngles,linspace(0,2*pi,nBins+1), ...
          'Normalization','pdf','FaceColor',"#9400D3", ...
          'FaceAlpha',0.3)
xlabel("azimuth angle [rad]")
ylabel("frequency density")
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels(["0" "\pi/2" "\pi" "3\pi/2" "2\pi"])
title("LOD azimuth angle")
axis tight
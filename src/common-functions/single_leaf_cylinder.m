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

function [Leaves,QSMbc] = single_leaf_cylinder(TargetDistributions, ...
                                               LeafProperties, ...
                                               leafArea, ...
                                               varargin)
% Default values
len = 0.5;
rad = 0.025;
inc = pi/4;
az  = -pi/4;
PetioleDirectionDistribution.flag = false;
Phyllotaxis.flag = false;

% Check additional parameters
i = 1;
NArg = numel(varargin);
while i <= NArg

    if ischar(varargin{i})

        switch lower(varargin{i})

            % Cylinder parameters
            case "cylinderproperties"
                CylProperties = varargin{i+1};
                len = CylProperties.cylinderLength;
                rad = CylProperties.cylinderRadius;
                inc = CylProperties.cylinderInclinationAngle;
                az  = CylProperties.cylinderAzimuthAngle;
                i = i + 1;

            % Petiole direction distribution
            case 'petioledirectiondistribution'
                assert(i < NArg && isa(varargin{i+1},'function_handle'),...
                       'Argument following ''PetioleDirectionDistribution'' should be a function handle.');
                PetioleDirectionDistribution.flag = true;
                PetioleDirectionDistribution.dist_fun = varargin{i+1};
                i = i + 1;

            % Phyllotaxis
             case 'phyllotaxis'
                assert(i < NArg && isa(varargin{i+1},'struct'), ...
                       'Argument following ''Phyllotaxis'' should be a struct.')
                Phyllotaxis = varargin{i+1};
                Phyllotaxis.flag = true;
                if PetioleDirectionDistribution.flag == true
                    warning('Petiole direction distribution cannot be used simultaneously with phyllotaxis enabled')
                end
                i = i + 1;

            otherwise
                warning(['Skipping unknown parameters: ''' varargin{i} '''']);

        end
    end
    i = i + 1;
end



% Generate QSM out of the cylinder
qsm.cylinder.start  = [0 0 0];
qsm.cylinder.axis   = (rotation_matrix([0 0 1],az)* ...
                       rotation_matrix([-1 0 0],inc)*[0 0 1]')';
qsm.cylinder.length = len;
qsm.cylinder.radius = rad;
qsm.cylinder.parent = 0;
qsm.cylinder.branch = 1;
QSMbc = QSMBCylindrical(qsm);

% Leaf area for the cylinder
ar = leafArea;

% Initialize leaf base area and candidate leaf area
LeavesInit = LeafModelTriangle(LeafProperties.vertices, ...
                               LeafProperties.triangles);

% Sample leaves from leaf size function
overSamplingFactor = 5;
leafScaleFactors = fun_leaf_size(overSamplingFactor*ar, ...
                                 LeavesInit.base_area, ...
                                 TargetDistributions.dTypeLSD,...
                                 TargetDistributions.fun_pLSD(0,0,0));
nLeaves = size(leafScaleFactors,1);
maxLeafSize = max(max(leafScaleFactors))*max(LeavesInit.base_dimensions);

% Attach the leaves to the cylinder with leaf orientation
% distribution
[leafDir,leafNormal,petioleStart,petioleEnd] = fun_leaf_orientation(...
    len,rad,inc,az, ...
    nLeaves, ...
    LeafProperties.petioleLengthLimits, ...
    TargetDistributions.dTypeLODinc, ...
    TargetDistributions.dTypeLODaz, ...
    TargetDistributions.fun_pLODinc(0,0,0), ...
    TargetDistributions.fun_pLODaz(0,0,0), ...
    PetioleDirectionDistribution, ...
    Phyllotaxis ...
    );

% Average leaf area
avgAr = mean(LeavesInit.base_area*(leafScaleFactors(:,1).^2));
% Estimate on total leaf count
totalLeafArea = ar;
leafCountEst = int64(1.1*round(totalLeafArea/avgAr));
% Initialize leaf object
Leaves = LeafModelTriangle(LeafProperties.vertices, ...
                           LeafProperties.triangles, ...
                           0,leafCountEst);

% Add leaves without intersections
Leaves = add_leaves_cylinder(Leaves,len,rad,inc,az,ar, ...
                             petioleStart,petioleEnd,leafDir,leafNormal,...
                             leafScaleFactors,maxLeafSize, ...
                             max(LeafProperties.petioleLengthLimits));

% Trim excess rows
Leaves.trim_slack;

end
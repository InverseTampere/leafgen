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

function Leaves = add_leaves_qsm(QSM,Leaves,leafScaleFactors,leafParent,...
                                 leafDir,leafNormal,petioleStart, ...
                                 petioleEnd,totalLeafArea)

%% Maximum leaf scaling and petiole length
maxLeafSize = max(max(leafScaleFactors))*max(Leaves.base_dimensions);
maxPetioleLen = max(sum((petioleEnd-petioleStart).^2,2));

%% Voxelization parameters

% Extreme values of cylinder QSM
treeBox = QSM.tree_limits;
minPoint = treeBox(1,:) - maxLeafSize - maxPetioleLen;
maxPoint = treeBox(2,:) + maxLeafSize + maxPetioleLen;

% Compute voxelization of cylinder QSM
QSMVoxelization = QSM.toVoxels(maxLeafSize, minPoint, maxPoint);

% Initialize voxelization of accepted leaves
LeafVoxelization = CubeVoxelization(maxLeafSize, minPoint, maxPoint);

%% Initialization of leaf transfromations

% Transform configuration matrix
TransformConfig = [...
%   x    y    z     x     y    z    x    y    z
    1    1    1     0     0  pi/4   0    0    0;
    1    1    1     0     0 -pi/4   0    0    0;
    1    1    1     0  pi/6    0    0    0    0;
    1    1    1     0 -pi/6    0    0    0    0;
    1    1    1  pi/6     0    0    0    0    0;
    1    1    1 -pi/6     0    0    0    0    0;
];

% Number of transforms.
nTransform = size(TransformConfig,1);

%% Initialization of necessary variables and parameters

% Number of candidate leaves
nLeafCandidate = size(leafDir,1);

% Number of configurations tried
nConfigsTried = nan(nLeafCandidate,1);

% Number of neighbours on initial run (1: cylinder, 2: leaf)
nNeighbour = nan(nLeafCandidate,2);

% Number of accepted leaves
nAccepted = 0;

% Total area of accepted leaves
areaAccepted = 0;

% Randomize the order of transform configurations
jTransform = randperm(nTransform);

% Flag for reaching target leaf area
targetReached = false;



%% Add leaves to the model

for iLeaf = 1:nLeafCandidate
    
    % Check if target area has been reached
    if areaAccepted >= totalLeafArea
        targetReached = true;
        break;
    end
    
    % Cube coordinates of leaf on last try
    lastCC = [];
    
    % Compute leaf neighbour at least once.
    computeLeafNeigbours = true;

    for iTransform = 1:nTransform+1
        
        % Paramenters of the current leaf
        origin = petioleEnd(iLeaf,:);
        dir    = leafDir(iLeaf,:);
        normal = leafNormal(iLeaf,:);
        scale  = leafScaleFactors(iLeaf,:);
        
        % Third axis (x) of leaf coordinate system
        side = cross(dir,normal);
        
        % Flag: neighbour index computation required
        computeNeighbours = false;
        
        % Do transformation if necessary
        if iTransform > 1
            
            % Transformation configuration
            config = TransformConfig(jTransform(iTransform-1),:);
            
            % Scaling: multiply current dimensions
            if any(config(1:3) ~= ones(1,3))
                scale = scale.*config(1:3);
            end
            
            % Elevation / around side axis
            if config(4)
                Rx = rotation_matrix(side,config(4));
                dir = (Rx*dir')';
                normal = (Rx*normal')';
            end
            
            % Rotation / around direction
            if config(5)
                Ry = rotation_matrix(dir,config(5));
                %side = (Ry*side')';
                normal = (Ry*normal')';
            end
            
            % Azimuth / around normal
            if config(6)
                Rz = rotation_matrix(normal,config(6));
                dir = (Rz*dir')';
                %side = (Rz*side')';
            end
            
            % Translation
            if any(config(7:9) ~= zeros(1,3))
                origin = origin + config(7:9);
            end
        end
        
        % Compute leaf center point
        cen = origin + 0.5*scale(2)*Leaves.base_dimensions(2)*dir;
        
        % Compute leaf triangles for intersection computations
        leafTris = Leaves.triangles(origin, dir, normal, scale);
        
        % Convert center point to cube voxelization coordinates
        leafCC = LeafVoxelization.get_coordinates(cen);
    
        % On first run neighbour computations are always required
        if iTransform == 1
            computeNeighbours = true;
            lastCC = leafCC; 
        % Otherwise check if cube coordinate of leaf has changed
        elseif any(leafCC ~= lastCC)
            computeNeighbours = true;
            lastCC = leafCC;
        end
        
        % Find neighbour blocks and their count only if changes have
        % occured
        if computeNeighbours
            % Get neighbour cylinder indices (exclude parent cylinder)
            cylNeighbour = QSMVoxelization.get_neighbor_objects(leafCC);

            % Number of cylinder neighbours
            nCylNei = length(cylNeighbour);
            
            if isnan(nNeighbour(iLeaf,1))
                nNeighbour(iLeaf,1) = nCylNei;
            end
        end

        % If any overlap, try next configuration
        if nCylNei && QSM.block_triangle_intersection(cylNeighbour,...
                                                        leafTris)
            continue;
        end
        
        % Find neighbour blocks and their count only if changes have
        % occured
        if computeNeighbours || computeLeafNeigbours
            % Get neighbour cylinder indices (exclude parent cylinder)
            leafNeighbour = LeafVoxelization.get_neighbor_objects(leafCC);
            
            % After first compute, respect <computeNeighbours>
            computeLeafNeigbours = false;

            % Number of cylinder neighbours
            nLeafNei = length(leafNeighbour);
            
            if isnan(nNeighbour(iLeaf,2))
                nNeighbour(iLeaf,2) = nLeafNei;
            end
        end

        % If any overlap, try next configuration
        if nLeafNei && Leaves.leaf_intersect(leafNeighbour,leafTris)
            continue;
        end
    
        % Set the parent of accepted leaf as NaN
        parent = leafParent(iLeaf);
        
        % Petiole start point
        petiole = petioleStart(iLeaf,:);

        % If inclination of normal is above pi/2, mirror normal to other
        % side
        if acos(dot(normal,[0,0,1])) > pi/2
            normal = -normal;
        end
        
        % No intersections: leaf is accepted and inserted into model
        [leafIndex,leafArea] = Leaves.add_leaf(origin,...
                                               dir,...
                                               normal,...
                                               scale,...
                                               parent,...
                                               petiole,...
                                               leafTris);
        
        % Add leaf to voxelization
        LeafVoxelization.add_object_by_cc(leafCC, leafIndex);
        
        % Increase count and total area.
        nAccepted = nAccepted + 1;
        areaAccepted = areaAccepted + leafArea;
        
        % Skip remaining transform configurations.
        break;
        
    end

    % Store applied transformation config id.
    nConfigsTried(iLeaf) = (iTransform - 1);
    
end

if targetReached
    iLeaf = iLeaf - 1;
end

% Clear possible empty rows.
Leaves.trim_slack();

end
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

function [Leaves,QSM] = transform_leaf_cylinders(qsm, ...
                                                 LeafCylinderLibrary, ...
                                                 TargetLADD, ...
                                                 ParamFunctions, ...
                                                 targetLeafArea)

%% Initialize struct-based QSM as QSMBCylindrical-class object

if isa(qsm,'struct')
    QSM = QSMBCylindrical(qsm);
elseif isa(qsm,'QSMBCylindrical')
    QSM = qsm;
else
    error("Invalid QSM type. Input QSM should be a struct or " + ...
          "a QSMBCylindrical-class object.")
end

%% Check the correctness of inputs
QSM = check_inputs_transform(QSM,TargetLADD,ParamFunctions,targetLeafArea); 
                             

%% Reading leaf cylinder library nodes

Nodes = LeafCylinderLibrary.Nodes;

%% Extract QSM properties

ParLabels = {'relative_height',...
             'relative_position',...
             'branch_order',...
             'axis',...
             'length',...
             'radius',...
             'is_last',...
             'start_point',...
             'branch_index',...
             'index_in_branch'};

CylinderParameters = QSM.get_block_properties(ParLabels);

% Preprocess the cylinder lengths to better suit the leaf cylinder library
lengthLimits = [Nodes.cylinderLength(1) Nodes.cylinderLength(end)];
[CylinderParameters,originalCylIndex] = preprocess_cylinders( ...
                                            CylinderParameters, ...
                                            lengthLimits, ...
                                            "absolute");

% Calculate the relative distance of cylinders along subbranch
CylinderParameters.relative_distance_along_subbranch = ...
    calculate_relative_distance_along_subbranch(CylinderParameters);

% Calculate compass direction of cylinders with respect to mean value of
% horizontal cylinder location
CylinderParameters.compass_direction = ...
    calculate_compass_direction(CylinderParameters);

%% Determining leaf area budgets for cylinders using LADD

% Relative leaf area budgets
relativeCylinderLeafArea = fun_leaf_area_density(CylinderParameters, ...
                                                 TargetLADD);

% Scale by target area to get cylinder leaf area budgets
cylinderLeafArea = targetLeafArea*relativeCylinderLeafArea;

%% Adding leaves to cylinders from the leaf-cylinder library

% Initialize LeafModelTriangle class object
Leaves = LeafModelTriangle(LeafCylinderLibrary.LeafBaseModel.vertices, ...
                           LeafCylinderLibrary.LeafBaseModel.triangles);

for iCyl = 1:length(cylinderLeafArea)

    % If cylinder leaf area budget is under 1cm^2 omit cylinder
    if cylinderLeafArea(iCyl) < 0.01^2
        continue
    end

    % LENGTH
    % Check if cylinder length exceeds or equals the max length node of 
    % library
    if CylinderParameters.length(iCyl) >= Nodes.cylinderLength(end)
        iLen = length(Nodes.cylinderLength);
    % Check if cylinder length is less than min length node of library
    elseif CylinderParameters.length(iCyl) < Nodes.cylinderLength(1)
        iLen = 1;
    % Find first length node below the value of cylinder length
    else 
        iLen = find(Nodes.cylinderLength > CylinderParameters.length(iCyl), ...
                    1,'first') - 1;
    end
    % Length scaling factor
    sfLen = CylinderParameters.length(iCyl)/Nodes.cylinderLength(iLen);


    % RADIUS
    % Check if cylinder radius exceeds or equals the max radius node of 
    % library
    if CylinderParameters.radius(iCyl) >= Nodes.cylinderRadius(end)
        iRad = length(Nodes.cylinderRadius);
    % Check if cylinder radius is less than min radius node of library
    elseif CylinderParameters.radius(iCyl) < Nodes.cylinderRadius(1)
        iRad = 1;
    % Find first radius node below the value of cylinder radius
    else 
        iRad = find(Nodes.cylinderRadius > CylinderParameters.radius(iCyl), ...
                    1,'first') - 1;
    end
    % Radius scaling factor
    sfRad = CylinderParameters.radius(iCyl)/Nodes.cylinderRadius(iRad);
    

    % INCLINATION ANGLE
    % Cylinder axis vector (unit length)
    cylAx = CylinderParameters.axis(iCyl,:) ...
            /norm(CylinderParameters.axis(iCyl,:));
    % Cylinder inclination angle
    cylIncAngle = acos(dot(cylAx,[0 0 1]));
    % Find the closest node to cylinder inc angle value
    [incNode,iInc] = find_closest_node(Nodes.cylinderInclinationAngle, ...
                                            cylIncAngle);
    % Inclination angle translation term
    trInc = cylIncAngle - incNode;


    % AZIMUTH ANGLE
    % Check if cylinder axis is parallel to vertical direction
    if norm(cylAx - [0 0 1]) < 1e-6 || norm(cylAx - [0 0 -1]) < 1e-6
        azNode = Nodes.cylinderAzimuthAngle(1);
        iAz = 1;
        trAz = 0;
    else
        % Cylinder axis projection on horizontal plane (unit length)
        cylAxHorz = [cylAx(1),cylAx(2),0]/norm([cylAx(1),cylAx(2),0]);
        % Cylinder azimuth angle
        if cylAx(1) > 0
            cylAzAngle = 2*pi - acos(dot(cylAxHorz,[0 1 0]));
        else
            cylAzAngle = acos(dot(cylAxHorz,[0 1 0]));
        end
        % Check if cylinder azimuth angle is between the first and last
        % node
        if cylAzAngle > Nodes.cylinderAzimuthAngle(end)
            if (cylAzAngle-Nodes.cylinderAzimuthAngle(end)) < (2*pi-cylAzAngle)
                azNode = Nodes.cylinderAzimuthAngle(end);
                iAz = length(Nodes.cylinderAzimuthAngle);
            else
                azNode = Nodes.cylinderAzimuthAngle(1);
                iAz = 1;
            end
        else
            % Find the closest node to cylinder az angle value
            [azNode,iAz] = find_closest_node(Nodes.cylinderAzimuthAngle, ...
                                             cylAzAngle);
        end
        % Azimuth angle translation term
        trAz = cylAzAngle - azNode;
    end


    % LEAF AREA
    % Check if cylinder leaf area exceeds the max leaf area node of library
    if cylinderLeafArea(iCyl) > Nodes.cylinderLeafArea(end)
        iAr = length(Nodes.cylinderLeafArea);
    else % Find first cylinder leaf area node above the value of cyl leaf area
        iAr = find(Nodes.cylinderLeafArea > cylinderLeafArea(iCyl), ...
                   1,'first');
    end
                                     

    % LEAF ORIENTATION DISTRIBUTION 
    % Leaf inclination angle distribution parameter values for the cylinder
    pInc = ParamFunctions.fun_pLODinc( ...
                CylinderParameters.relative_height, ...
                CylinderParameters.relative_distance_along_subbranch, ...
                CylinderParameters.compass_direction);
    % Add dummy element for constant inclination angle
    if isscalar(pInc)
        pInc = [pInc(1) 0];
    end
    % Leaf azimuth angle distribution parameter values for the cylinder
    pAz  = ParamFunctions.fun_pLODaz( ...
                CylinderParameters.relative_height, ...
                CylinderParameters.relative_distance_along_subbranch, ...
                CylinderParameters.compass_direction);
    % Add dummy element for constant azimuth angle
    if isscalar(pAz)
        pAz = [pAz(1) 0];
    end
    % LOD nodes
    [~,iLodInc1] = find_closest_node(Nodes.pLODinc1,pInc(1));
    [~,iLodInc2] = find_closest_node(Nodes.pLODinc2,pInc(2));
    [~,iLodAz1]  = find_closest_node(Nodes.pLODaz1,pAz(1));
    [~,iLodAz2]  = find_closest_node(Nodes.pLODaz2,pAz(2));

    % LEAF SIZE DISTRIBUTION
    % Leaf size distribution parameter values for the cylinder
    pSize = ParamFunctions.fun_pLSD( ...
                CylinderParameters.relative_height, ...
                CylinderParameters.relative_distance_along_subbranch, ...
                CylinderParameters.compass_direction);
    % Add dummy element for constant leaf size
    if isscalar(pSize)
        pSize = [pSize(1) 0];
    end
    % LSD nodes
    [~,iLsd1] = find_closest_node(Nodes.pLSD1,pSize(1));
    [~,iLsd2] = find_closest_node(Nodes.pLSD2,pSize(2));


    % PICK CORRESPONDING LEAF OBJECT FROM THE LIBRARY
    iNodeObj = randi(LeafCylinderLibrary.Properties.nLeafObjectsPerNode);
    iNode = sub2ind(LeafCylinderLibrary.Properties.nNodesPerLibVar, ...
                    iLen,iRad,iInc,iAz,iAr,iLodInc1,iLodInc2,iLodAz1, ...
                    iLodAz2,iLsd1,iLsd2,iNodeObj);
    LibraryObj = LeafCylinderLibrary.LeafObjects(iNode).Leaves;

    % TRANSFORM LEAF PARAMETERS TO CORRESPOND THE QSM CYLINDER SIZE AND
    % ORIENTATION
    % Library object axis
    if incNode == 0
        libObjAxis = [0 0 1];
    elseif incNode == pi
        libObjAxis = [0 0 -1];
    else
        libObjHorzAxis = [cos(azNode+pi/2) sin(azNode+pi/2)];
        libObjHorzAxis = libObjHorzAxis/norm(libObjHorzAxis);
        libObjAxis = [abs(sin(incNode))*libObjHorzAxis cos(incNode)];
        libObjAxis = libObjAxis/norm(libObjAxis);
    end
    % Vector from petiole start to leaf start
    petioleToLeaf = LibraryObj.leaf_start_point ...
                    - LibraryObj.petiole_start_point;
    % Petiole start point axis-wise translation
    petioleStartPoints = (sfLen-1)* ...
                         (LibraryObj.petiole_start_point*libObjAxis') ...
                         .*libObjAxis + LibraryObj.petiole_start_point;
    % Petiole start point radius-wise translation
    radVecPetiole = petioleStartPoints - (petioleStartPoints*libObjAxis') ...
                 .*libObjAxis;
    petioleStartPoints = (petioleStartPoints*libObjAxis').*libObjAxis ...
                      + sfRad*radVecPetiole;
    % Translated leaf start points
    leafStartPoints = petioleStartPoints + petioleToLeaf;

    % Leaf and petiole start point and leaf direction and normal rotation 
    % around the origin
    leafDirections = LibraryObj.leaf_direction;
    leafNormals = LibraryObj.leaf_normal;
    if norm(libObjAxis - [0 0 1] ) > 1e-6 && ...
       norm(libObjAxis - [0 0 -1]) > 1e-6
        % Inclination angle rotation matrix
        rotAxis = cross([0 0 1],libObjAxis);
        rotAxis = rotAxis/norm(rotAxis);
        incRM = rotation_matrix(rotAxis,trInc);
        % Azimuth angle rotation matrix
        azRM = rotation_matrix([0 0 1],trAz);
    else
        % Inclination angle rotation matrix
        rotAxis = [-1 0 0];
        incRM = rotation_matrix(rotAxis,trInc);
        % Azimuth angle rotation matrix (no azimuth rotation
        azRM = eye(3);
    end    
    % Inclination angle rotations
    leafStartPoints = (incRM*leafStartPoints')';
    petioleStartPoints = (incRM*petioleStartPoints')';
    leafDirections  = (incRM*leafDirections')';
    leafNormals     = (incRM*leafNormals')';
    % Azimuth angle rotations
    leafStartPoints = (azRM*leafStartPoints')';
    petioleStartPoints = (azRM*petioleStartPoints')';
    leafDirections  = (azRM*leafDirections')';
    leafNormals     = (azRM*leafNormals')';

    % Fix leaf inclination angles of over pi/2 radians possibly caused by
    % the rotations
    iIncFix = find(leafNormals(:,3) < 0);
    leafNormals(iIncFix,:) = -leafNormals(iIncFix,:);
    
    % Leaf parent cylinder
    leafParent = originalCylIndex(iCyl);

    % Library leaf object base area and leaf area scalings
    leafBaseArea  = LibraryObj.base_area;
    leafScales    = LibraryObj.leaf_scale;

    % Check if cylinder leaf area budget exceeds the leaf area of
    % library object
    if cylinderLeafArea(iCyl) > LibraryObj.leaf_area
        endIndex = LibraryObj.leaf_count;
    elseif  cylinderLeafArea(iCyl) > 0.5*leafBaseArea*(leafScales(1,1)^2)
        % Reduce leaves to correspond to the leaf area budget of the cyl
        remainingArea = LibraryObj.leaf_area;
        nLeavesInObj = LibraryObj.leaf_count;
        k = 0;
        nExcessLeaves = 0;
        while remainingArea > cylinderLeafArea(iCyl) && k < nLeavesInObj
            remainingArea = remainingArea - leafBaseArea ...
                            *(leafScales(end-k,1)^2);
            k = k + 1;
            nExcessLeaves = k - 1;
        end
        endIndex = LibraryObj.leaf_count - nExcessLeaves;
    else % leaf area budget is under 0.5 times leaf size -> omit cylinder
        continue
    end

    % TRANSLATE THE LEAVES TO THE CORRECT POSITION
    leafStartPoints = CylinderParameters.start_point(iCyl,:) ...
                      + leafStartPoints;
    petioleStartPoints = CylinderParameters.start_point(iCyl,:) ...
                         + petioleStartPoints;

    % ADD LEAVES TO THE TREE-LEAF MODEL
    if iCyl == 1 && endIndex > 0
        Leaves.leaf_count = endIndex;
        Leaves.leaf_area = sum(leafBaseArea*(leafScales(1:endIndex,1).^2));
        Leaves.leaf_start_point    = leafStartPoints(1:endIndex,:);
        Leaves.leaf_scale          = leafScales(1:endIndex,:);
        Leaves.leaf_direction      = leafDirections(1:endIndex,:);
        Leaves.leaf_normal         = leafNormals(1:endIndex,:);
        Leaves.leaf_parent         = leafParent*ones(endIndex,1);
        Leaves.petiole_start_point = petioleStartPoints(1:endIndex,:);
    elseif endIndex > 0
        Leaves.leaf_count = Leaves.leaf_count + endIndex;
        Leaves.leaf_area           = Leaves.leaf_area ...
                                     + sum(leafBaseArea ...
                                           *(leafScales(1:endIndex,1).^2));
        Leaves.leaf_start_point    = cat(1,Leaves.leaf_start_point, ...
                                      leafStartPoints(1:endIndex,:));
        Leaves.leaf_scale          = cat(1,Leaves.leaf_scale, ...
                                      leafScales(1:endIndex,:));
        Leaves.leaf_direction      = cat(1,Leaves.leaf_direction, ...
                                      leafDirections(1:endIndex,:));
        Leaves.leaf_normal         = cat(1,Leaves.leaf_normal, ...
                                      leafNormals(1:endIndex,:));
        Leaves.leaf_parent         = cat(1,Leaves.leaf_parent, ...
                                      leafParent*ones(endIndex,1));
        Leaves.petiole_start_point = cat(1,Leaves.petiole_start_point, ...
                                      petioleStartPoints(1:endIndex,:));
    end

clearvars LibraryObj

end

% If target area is surpassed, remove leaves until target is reached
if Leaves.leaf_area > targetLeafArea
    % Randomize leaf order
    randOrder = randperm(Leaves.leaf_count);
    Leaves.leaf_start_point = Leaves.leaf_start_point(randOrder,:);
    Leaves.leaf_scale = Leaves.leaf_scale(randOrder,:);
    Leaves.leaf_direction = Leaves.leaf_direction(randOrder,:);
    Leaves.leaf_normal = Leaves.leaf_normal(randOrder,:);
    Leaves.leaf_parent = Leaves.leaf_parent(randOrder);
    Leaves.petiole_start_point = Leaves.petiole_start_point(randOrder,:);
    % Remove leaves until just above target area
    reducedArea = Leaves.leaf_area;
    k = 0;
    nRemove = 0;
    while reducedArea > targetLeafArea
        reducedArea = reducedArea - ...
                      Leaves.base_area*Leaves.leaf_scale(end-k,1).^2;
        k = k + 1;
        nRemove = k - 1;
    end
    newEndInd = Leaves.leaf_count - nRemove;
    Leaves.leaf_count = newEndInd;
    Leaves.leaf_area = Leaves.leaf_area - ...
                       Leaves.base_area* ...
                          sum(Leaves.leaf_scale((end-nRemove+1):end,1).^2);
    Leaves.leaf_start_point = Leaves.leaf_start_point(1:newEndInd,:);
    Leaves.leaf_scale = Leaves.leaf_scale(1:newEndInd,:);
    Leaves.leaf_direction = Leaves.leaf_direction(1:newEndInd,:);
    Leaves.leaf_normal = Leaves.leaf_normal(1:newEndInd,:);
    Leaves.leaf_parent = Leaves.leaf_parent(1:newEndInd);
    Leaves.petiole_start_point = Leaves.petiole_start_point(1:newEndInd,:);
end

% Display a warning if generated leaf area is lower than target leaf area
if Leaves.leaf_area/targetLeafArea < 0.99
    warning("Final leaf area " + num2str(Leaves.leaf_area) + " m^2 is " ...
            + num2str(targetLeafArea-Leaves.leaf_area) + " m^2 lower " ...
            + "than the target leaf area " + num2str(targetLeafArea) ...
            + " m^2. This is due to cylinder leaf area budgets " ...
            + "assigned by LADD surpassing the maximum leaf area of " ...
            + "the leaf cylinder library cylinders.")
end

end
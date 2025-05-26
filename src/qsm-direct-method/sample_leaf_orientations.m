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

function [leafDir,leafNor,petioleS,petioleE] = sample_leaf_orientations(...
                                              CylinderParameters, ...
                                              leafParent, ...
                                              TargetDistributions, ...
                                              LeafProperties, ...
                                              PetiolegDirectionDistribution,...
                                              Phyllotaxis)

% Petiole length limits
petioleLengthLimits = LeafProperties.petioleLengthLimits;

nLeavesTotal = size(leafParent,1);
leafDir = zeros(nLeavesTotal,3);
leafNor = zeros(nLeavesTotal,3);
petioleS   = zeros(nLeavesTotal,3);
petioleE   = zeros(nLeavesTotal,3);
iLeaf = 1;
for iCyl = 1:length(CylinderParameters.relative_height)
    % Number of leaves on cylinder
    nLeavesCyl = sum(leafParent == iCyl);
    if nLeavesCyl > 0
        % Cylinder structural variable coordinates
        svCoords = [CylinderParameters.relative_height(iCyl), ...
             CylinderParameters.relative_distance_along_subbranch(iCyl),...
             CylinderParameters.compass_direction(iCyl)];
        % Cylinder length
        len = CylinderParameters.length(iCyl);
        % Cylinder radius
        rad = CylinderParameters.radius(iCyl);
        % Cylinder axis
        axis = CylinderParameters.axis(iCyl,:) ...
               /norm(CylinderParameters.axis(iCyl,:));
        % Cylinder inclination angle
        inc = acos(dot(axis,[0 0 1]));
        % Cylinder azimuth angle (counterclockwise wrt. north direction)
        if dot(axis,[0 0 1]) > 0.999 % vertical axis
            az = 0;
        else
            axisXY = [axis(1) axis(2) 0];
            axisXY = axisXY/norm(axisXY);
            if axis(1) <= 0
                az = acos(dot(axisXY,[0 1 0]));
            else
                az = 2*pi - acos(dot(axisXY,[0 1 0]));
            end
        end
    
        % LOD distribution types and parameters
        dTypeLODinc = TargetDistributions.dTypeLODinc;
        dTypeLODaz  = TargetDistributions.dTypeLODaz;
        dParametersLODinc = TargetDistributions.fun_pLODinc(svCoords(1),...
                                                            svCoords(2),...
                                                            svCoords(3));
        dParametersLODaz  = TargetDistributions.fun_pLODaz(svCoords(1),...
                                                           svCoords(2),...
                                                           svCoords(3));
        % Sample leaf orientation and petiole coordinates
        [dir,nor,ps,pe] = fun_leaf_orientation(len,rad,inc,az, ...
                                              nLeavesCyl, ...
                                              petioleLengthLimits, ...
                                              dTypeLODinc, ...
                                              dTypeLODaz, ...
                                              dParametersLODinc, ...
                                              dParametersLODaz, ...
                                              PetiolegDirectionDistribution,...
                                              Phyllotaxis);
        nLeavesAdded = size(dir,1);
        leafDir(iLeaf:(iLeaf+(nLeavesAdded-1)),:) = dir;
        leafNor(iLeaf:(iLeaf+(nLeavesAdded-1)),:) = nor;
        petioleS(iLeaf:(iLeaf+(nLeavesAdded-1)),:) = ...
                              ps + CylinderParameters.start_point(iCyl,:);
        petioleE(iLeaf:(iLeaf+(nLeavesAdded-1)),:) = ...
                              pe + CylinderParameters.start_point(iCyl,:);

        iLeaf = iLeaf + nLeavesAdded;
    end
end

end
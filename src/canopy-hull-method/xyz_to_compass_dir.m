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

function compassDir = xyz_to_compass_dir(leafStartPoint,stemCoordinates)
% Index of the first stem coordinate with z-value above the leaf height
if stemCoordinates(end,3) < leafStartPoint(3)
    iSC = size(stemCoordinates,1);
else
    iSC = find(stemCoordinates(:,3) > leafStartPoint(3),1);
end
% Relative z-position of leaf between the stem coordinate nodes
relPos = (leafStartPoint(3)- stemCoordinates(iSC-1,3)) ...
          /(stemCoordinates(iSC,3) - stemCoordinates(iSC-1,3));
% Location of the stem center on the height of the leaf
stemCen = relPos*(stemCoordinates(iSC,:) - stemCoordinates(iSC-1,:)) ...
          + stemCoordinates(iSC-1,:);
% Horizontal distance from stem to leaf
stemToLeaf = [leafStartPoint(1:2) 0] - [stemCen(1:2) 0];
distStemToLeaf = norm(stemToLeaf);
% Counter-clockwise angle between north and leaf direction wrt. the stem
unitStemToLeaf = stemToLeaf/distStemToLeaf;
if unitStemToLeaf(1) < 0
    compassDir = acos(dot(unitStemToLeaf,[0 1 0]));
else
    compassDir = 2*pi - acos(dot(unitStemToLeaf,[0 1 0]));
end
end
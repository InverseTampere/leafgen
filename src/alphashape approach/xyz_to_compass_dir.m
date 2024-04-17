function compassDir = xyz_to_compass_dir(leafStartPoint,stemCoordinates)
% Index of stem coordinate with z-value below the leaf height
iSC = find(stemCoordinates(:,3) > leafStartPoint(3),1);
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
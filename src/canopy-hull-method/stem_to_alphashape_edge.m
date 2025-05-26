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

function edgeDistance = stem_to_alphashape_edge(shp, ...
                                                heightValue, ...
                                                compassValue, ...
                                                maxHorzDist, ...
                                                maxHeight, ...
                                                stemCoordinates)

% Probing the edge of alpha shape
nPP = 100;
yCoord = 2*maxHorzDist*linspace(0,1,nPP)';
initPP = [zeros(nPP,1) yCoord zeros(nPP,1)];
probePoints = (rotation_matrix([0 0 1],compassValue)*initPP')' ...
    + [zeros(nPP,2) maxHeight*heightValue*ones(nPP,1)];
if stemCoordinates(end,3) < maxHeight*heightValue
    iSC = size(stemCoordinates,1);
else
    iSC = find(stemCoordinates(:,3) > maxHeight*heightValue,1);
end
relPos = (maxHeight*heightValue-stemCoordinates(iSC-1,3)) ...
    /(stemCoordinates(iSC,3)-stemCoordinates(iSC-1,3));
stemCen = relPos*(stemCoordinates(iSC,:) - stemCoordinates(iSC-1,:)) ...
    + stemCoordinates(iSC-1,:);
probePoints = probePoints + [stemCen(1:2) 0];
tf = inShape(shp,probePoints);
edgeIndex = find(~tf,1,'first') - 1;
if edgeIndex < 2
    % The first two probe points are not inside the shape
    edgeDistance = 0;
else
    edgeDistance = yCoord(edgeIndex);
end

end
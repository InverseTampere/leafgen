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
iSC = find(stemCoordinates(:,3) > maxHeight*heightValue,1);
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
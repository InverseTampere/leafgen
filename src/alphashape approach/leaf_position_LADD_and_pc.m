function leafStartPoint = leaf_position_LADD_and_pc(aShape, ...
                                                    fDist_h, ...
                                                    maxfDist_h, ...
                                                    xEdges, ...
                                                    yEdges, ...
                                                    zEdges, ...
                                                    pcVoxels, ...
                                                    voxelInd, ...
                                                    maxHeight ...
                                                    )
accepted = 0;
while accepted == 0
    % Proposal values for relative height
    hProposal = rand(1);
    % LADD value on propsal point
    funValue = fDist_h(hProposal);
    vertValue = rand(1)*maxfDist_h;
    if all(vertValue < funValue)
        hLeaf = maxHeight*hProposal;
        accepted = 1;
    end
end

% Corresponding height index in voxelization
iz = find(zEdges>hLeaf,1,'first') - 1;
accepted = 0;
zSlice = pcVoxels(:,:,iz);
zSliceInds = voxelInd(:,:,iz);
zSliceCol = zSlice(:);
zSliceIndsCol = zSliceInds(:);
cumulativeSum = cumsum(zSliceCol);
cumulativePDF = cumulativeSum./cumulativeSum(end);
while accepted == 0
    k = find(cumulativePDF>rand(1),1,'first');
    xyInds = zSliceIndsCol{k};
    ix = xyInds(1);
    iy = xyInds(2);
    % Uniformly sample a point inside the voxel
    proposal = rand(1,3).*[xEdges(ix+1)-xEdges(ix) ...
                           yEdges(iy+1)-yEdges(iy) ...
                           0] ...
               + [xEdges(ix) yEdges(iy) hLeaf];
    if inShape(aShape,proposal)
        leafStartPoint = proposal;
        accepted = 1;
    end
end
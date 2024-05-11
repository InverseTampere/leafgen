function leafStartPoint = sample_leaf_position(aShape, ...
                              fDist_h,fDist_d,fDist_c, ...
                              maxfDist_h,maxfDist_d,maxfDist_c, ...
                              maxHorzDist,maxHeight, ...
                              stemCoordinates, pcSamplingBins, ...
                              pcSamplingWeights, ...
                              xEdges,yEdges,zEdges, ...
                              pcVoxels,voxelInd,voxelCenDir)




% Use uniform sampling of alphashape (or point cloud) in positioning
if isempty(fDist_h) && isempty(fDist_d) && isempty(fDist_c)

    relHeight = rand(1);
    iPCS = find(relHeight<pcSamplingBins,1,'first') - 1;
    if pcSamplingWeights(iPCS) >= rand(1)
        % Point cloud sampling
        cumulativeSum = cumsum(pcVoxels(:));
        cumulativePDF = cumulativeSum./cumulativeSum(end);
        k = find(cumulativePDF>rand(1),1,'first');
        xyzInds = voxelInd{k};
        ix = xyzInds(1);
        iy = xyzInds(2);
        iz = xyzInds(3);
        % Uniformly sample a point inside the voxel
        proposal = rand(1,3).*[xEdges(ix+1)-xEdges(ix) ...
                               yEdges(iy+1)-yEdges(iy) ...
                               zEdges(iz+1)-zEdges(iz)] ...
                   + [xEdges(ix) yEdges(iy) zEdges(iz)];
        leafStartPoint = proposal;
    else
        % Uniform sampling inside alphashape
        proposal = [0 0 relHeight*maxHeight];
        accepted = 0;
        k = 0;
        while accepted == 0
            proposal(1) = rand(1)*(xEdges(end)-xEdges(1)) + xEdges(1);
            proposal(2) = rand(1)*(yEdges(end)-yEdges(1)) + yEdges(1);
            if inShape(aShape,proposal)
                leafStartPoint = proposal;
                accepted = 1;
            end
            k = k + 1;
            if k > 10000
                % No alphashape found on height level, resample height
                proposal(3) = rand(1)*maxHeight;
                k = 0;
            end
        end
        
    end



% Use heightwise LADD and uniform sampling of alphashape (or point cloud) 
% in positioning
elseif isempty(fDist_d) && isempty(fDist_c)

    accepted = 0;
    while accepted == 0
        % Sample height value from heightwise LADD
        hLeaf = rejection_sampling(fDist_h,maxfDist_h);
        % Find index of point cloud sampling height bin
        iPCS = find(hLeaf<pcSamplingBins,1,'first') - 1;
        if pcSamplingWeights(iPCS) >= rand(1)
            % Corresponding height index in voxelization
            iz = find(zEdges>hLeaf,1,'first') - 1;
            zSlice = pcVoxels(:,:,iz);
            zSliceInds = voxelInd(:,:,iz);
            zSliceCol = zSlice(:);
            zSliceIndsCol = zSliceInds(:);
            cumulativeSum = cumsum(zSliceCol);
            cumulativePDF = cumulativeSum./cumulativeSum(end);
            k = find(cumulativePDF>rand(1),1,'first');
            xyInds = zSliceIndsCol{k};
            ix = xyInds(1);
            iy = xyInds(2);
            % Uniformly sample a point inside the voxel
            proposal = rand(1,3).*[xEdges(ix+1)-xEdges(ix) ...
                                   yEdges(iy+1)-yEdges(iy) ...
                                   0] ...
                       + [xEdges(ix) yEdges(iy) hLeaf];
            leafStartPoint = proposal;
        else
            % Uniform sampling inside alphashape
            proposal = [0 0 hLeaf*maxHeight];
            k = 0;
            while accepted == 0
                proposal(1) = rand(1)*(xEdges(end)-xEdges(1)) + xEdges(1);
                proposal(2) = rand(1)*(yEdges(end)-yEdges(1)) + yEdges(1);
                if inShape(aShape,proposal)
                    leafStartPoint = proposal;
                    accepted = 1;
                end
                k = k + 1;
                if k > 10000
                    % No alphashape found on height level, resample height
                    break
                end
            end
        end
    end
    


% Use heightwise and directionwise LADD (or point cloud) in positioning
elseif isempty(fDist_d)


    accepted = 0;
    while accepted == 0
        % Proposal values for relative height and compass direction
        hProposal = rand(1);
        cProposal = 2*pi*rand(1);
        % LADD values on propsal point
        funValues = [fDist_h(hProposal) fDist_c(cProposal)];
        vertValues = rand(1,2).*[maxfDist_h maxfDist_c];
        if any(vertValues > funValues)
            continue
        end
        hLeaf = maxHeight*hProposal;
        cLeaf = cProposal;

        % Find index of point cloud sampling height bin
        iPCS = find(hProposal<pcSamplingBins,1,'first') - 1;
        if pcSamplingWeights(iPCS) >= rand(1) 
            % Corresponding height index in voxelization
            iz = find(zEdges>hLeaf,1,'first') - 1;
            zSlice = pcVoxels(:,:,iz);
            zSliceCol = zSlice(:);
            zSliceInds = voxelInd(:,:,iz);
            zSliceIndsCol = zSliceInds(:);
            zSliceCenDir = voxelCenDir(:,:,iz);
            zSliceCenDirCol = zSliceCenDir(:);
            % Corresponding voxels around the compass direction
            sectorWidth = 2*pi/10;
            inSector = all([(zSliceCenDirCol < cLeaf+0.5*sectorWidth) ...
                (zSliceCenDirCol > cLeaf-0.5*sectorWidth)], ...
                2);
            sectorSliceIndsCol = zSliceIndsCol(inSector);
            sectorSliceCol = zSliceCol(inSector);
            if sum(sectorSliceCol) == 0
                continue
            end
            cumulativeSum = cumsum(sectorSliceCol);
            cumulativePDF = cumulativeSum./cumulativeSum(end);
            k = find(cumulativePDF>rand(1),1,'first');
            xyInds = sectorSliceIndsCol{k};
            ix = xyInds(1);
            iy = xyInds(2);
            % Uniformly sample a point inside the voxel
            proposal = rand(1,3).*[xEdges(ix+1)-xEdges(ix) ...
                                   yEdges(iy+1)-yEdges(iy) ...
                                   0] ...
                       + [xEdges(ix) yEdges(iy) hLeaf];
            leafStartPoint = proposal;
            accepted = 1;
        else
            % Probing the edge of alpha shape
            edgeValue = stem_to_alphashape_edge(aShape, ...
                                                hProposal, ...
                                                cProposal, ...
                                                maxHorzDist, ...
                                                maxHeight, ...
                                                stemCoordinates);
            % Uniformly sample a point inside the alphashape in the
            % defined direction
            dLeaf = rand(1)*edgeValue;
            % Stem center on the corresponding height
            iSC = find(stemCoordinates(:,3) > hLeaf,1);
            relPos = (hLeaf-stemCoordinates(iSC-1,3)) ...
                     /(stemCoordinates(iSC,3)-stemCoordinates(iSC-1,3));
            stemCen = relPos*(stemCoordinates(iSC,:) ...
                              - stemCoordinates(iSC-1,:)) ...
                      + stemCoordinates(iSC-1,:);
            % Calculate the leaf start point
            leafStartPoint = (rotation_matrix([0 0 1],cLeaf) ...
                              *[0 1 0]')'*dLeaf ...
                             +[stemCen(1:2) hLeaf];
            accepted = 1;
        end
    end



% Use only LADD in positioning
elseif isa(fDist_h,"function_handle") && isa(fDist_d,"function_handle") ...
       && isa(fDist_c,"function_handle")

    accepted = 0;
    while accepted == 0
        % Proposal values for the variables
        hProposal = rand(1);
        dProposal = rand(1);
        cProposal = 2*pi*rand(1);
        % LADD value on propsal point
        funValues = [fDist_h(hProposal),fDist_d(dProposal), ...
                     fDist_c(cProposal)];
        vertValues = rand(1,3).*[maxfDist_h,maxfDist_d,maxfDist_c];
        if all(vertValues < funValues)
            % Probing the edge of alpha shape
            edgeValue = stem_to_alphashape_edge(aShape, ...
                                                hProposal, ...
                                                cProposal, ...
                                                maxHorzDist, ...
                                                maxHeight, ...
                                                stemCoordinates);
            % Check the if the distance to edge is nonzero
            if edgeValue == 0
                % The first two probe points are not inside the shape,
                % sample a new proposal point
                continue
            end
            % Position coordinates of the leaf start point
            hLeaf = maxHeight*hProposal;
            dLeaf = edgeValue*dProposal;
            cLeaf = cProposal;
            % Stem center on the corresponding height
            iSC = find(stemCoordinates(:,3) > hLeaf,1);
            relPos = (hLeaf-stemCoordinates(iSC-1,3)) ...
                     /(stemCoordinates(iSC,3)-stemCoordinates(iSC-1,3));
            stemCen = relPos*(stemCoordinates(iSC,:) ...
                              - stemCoordinates(iSC-1,:)) ...
                      + stemCoordinates(iSC-1,:);
            % Calculate the leaf start point
            leafStartPoint = (rotation_matrix([0 0 1],cLeaf) ...
                              *[0 1 0]')'*dLeaf ...
                             +[stemCen(1:2) hLeaf];
            accepted = 1;
        end
    end
    
end
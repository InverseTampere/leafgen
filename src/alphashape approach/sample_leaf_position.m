function [leafSP,PCSampling] = sample_leaf_position(aShape, ...
                                   fDist_h,fDist_d,fDist_c, ...
                                   maxfDist_h,maxfDist_d,maxfDist_c, ...
                                   xMin,xMax,yMin,yMax,maxHeight, ...
                                   stemCoordinates,PCSampling)


% Use uniform sampling of alphashape (or point cloud) in positioning
if isempty(fDist_h) && isempty(fDist_d) && isempty(fDist_c)

    uniformSampling = true;
    relHeight = rand(1);
    iPCS = find(relHeight<PCSampling.binEdges,1,'first') - 1;
    if PCSampling.ratio(iPCS) <= PCSampling.weights(iPCS) ...
            && PCSampling.flag == true
        % Voxels inside the point cloud sampling bin
        izLB = find(PCSampling.zEdges>maxHeight ...
                    *PCSampling.binEdges(iPCS),1,'first') - 1;
        izUB = find(PCSampling.zEdges>=maxHeight ...
                    *PCSampling.binEdges(iPCS+1),1,'first') - 1;
        % Cumulative sum of voxels inside sampling bin
        samplingBinVoxels = PCSampling.voxels(:,:,izLB:izUB);
        cumulativeSum = cumsum(samplingBinVoxels(:));
        samplingBinVoxelInds = PCSampling.voxelIndex(:,:,izLB:izUB);
        samplingBinVoxelIndsCol = samplingBinVoxelInds(:);
        % Check whether there are point cloud voxels inside the sampling
        % bin
        if cumulativeSum(end) > 0 
            uniformSampling = false;
            % Point cloud sampling
            cumulativePDF = cumulativeSum./cumulativeSum(end);
            k = find(cumulativePDF>rand(1),1,'first');
            % if izLB > 1
            %     belowVoxels = PCSampling.voxels(:,:,1:(izLB-1));
            %     k = k + length(belowVoxels(:));
            % end
            % xyzInds = PCSampling.voxelIndex{k};
            xyzInds = samplingBinVoxelIndsCol{k};
            ix = xyzInds(1);
            iy = xyzInds(2);
            iz = xyzInds(3);
            % Uniformly sample a point inside the voxel
            proposal = rand(1,3) ...
                       .*[PCSampling.xEdges(ix+1)-PCSampling.xEdges(ix) ...
                          PCSampling.yEdges(iy+1)-PCSampling.yEdges(iy) ...
                          PCSampling.zEdges(iz+1)-PCSampling.zEdges(iz)]...
                       + [PCSampling.xEdges(ix) PCSampling.yEdges(iy) ...
                          PCSampling.zEdges(iz)];
            % Set the sampled point as leaf start point
            leafSP = proposal;
            % Increment sampling counters and calculate new ratio
            PCSampling.nPCSampled(iPCS) = PCSampling.nPCSampled(iPCS) + 1;
            PCSampling.nTotalSampled(iPCS) = ...
                                       PCSampling.nTotalSampled(iPCS) + 1;
            PCSampling.ratio(iPCS) = PCSampling.nPCSampled(iPCS) ...
                                     /PCSampling.nTotalSampled(iPCS);
        end
    end

    if uniformSampling == true
        proposal = [0 0 0];
        % Sample the position uniformly
        accepted = 0;
        while accepted == 0
            proposal(1) = rand(1)*(xMax-xMin) + xMin;
            proposal(2) = rand(1)*(yMax-yMin) + yMin;
            proposal(3) = rand(1)*maxHeight;
            % Check if proposed point is inside the alpha shape
            if inShape(aShape,proposal)
                % Set the sampled point as leaf start point
                leafSP = proposal;
                % Increment sampling counter and calculate new ratio
                rh = proposal(3)/maxHeight;
                iPCS = find(rh<PCSampling.binEdges,1,'first') - 1;
                PCSampling.nTotalSampled(iPCS) = ...
                                       PCSampling.nTotalSampled(iPCS) + 1;
                PCSampling.ratio(iPCS) = PCSampling.nPCSampled(iPCS) ...
                                         /PCSampling.nTotalSampled(iPCS);
                accepted = 1;
            end
        end
        
    end



% Use heightwise LADD and uniform sampling of alphashape (or point cloud) 
% in positioning
elseif isempty(fDist_d) && isempty(fDist_c)
    
    uniformSampling = false;
    accepted = 0;
    while accepted == 0
        % Sample height value from heightwise LADD
        hLeaf = rejection_sampling(fDist_h,maxfDist_h);
        % Find index of point cloud sampling height bin
        iPCS = find(hLeaf<PCSampling.binEdges,1,'first') - 1;
        if all([PCSampling.ratio(iPCS) <= PCSampling.weights(iPCS), ...
                uniformSampling == false, ...
                PCSampling.flag == true])
            % Corresponding height index in voxelization
            iz = find(PCSampling.zEdges>hLeaf,1,'first') - 1;
            % Find all voxels for the height level
            zSlice = PCSampling.voxels(:,:,iz);
            zSliceInds = PCSampling.voxelIndex(:,:,iz);
            zSliceCol = zSlice(:);
            zSliceIndsCol = zSliceInds(:);
            % If no voxels on height level, do uniform sampling
            if sum(zSliceCol) == 0
                uniformSampling = true;
                continue
            end
            cumulativeSum = cumsum(zSliceCol);
            cumulativePDF = cumulativeSum./cumulativeSum(end);
            k = find(cumulativePDF>rand(1),1,'first');
            xyInds = zSliceIndsCol{k};
            ix = xyInds(1);
            iy = xyInds(2);
            % Uniformly sample a point inside the voxel
            proposal = rand(1,3) ...
                       .*[PCSampling.xEdges(ix+1)-PCSampling.xEdges(ix) ...
                          PCSampling.yEdges(iy+1)-PCSampling.yEdges(iy) ...
                          0] ...
                       + [PCSampling.xEdges(ix) PCSampling.yEdges(iy) ...
                          hLeaf];
            % Set the sampled point as leaf start point
            leafSP = proposal;
            % Increment sampling counters and calculate new ratio
            PCSampling.nPCSampled(iPCS) = PCSampling.nPCSampled(iPCS) + 1;
            PCSampling.nTotalSampled(iPCS) = ...
                                       PCSampling.nTotalSampled(iPCS) + 1;
            PCSampling.ratio(iPCS) = PCSampling.nPCSampled(iPCS) ...
                                     /PCSampling.nTotalSampled(iPCS);
            accepted = 1;
        end

        if uniformSampling == true
            % Uniform sampling inside alphashape
            proposal = [0 0 hLeaf*maxHeight];
            while accepted == 0
                proposal(1) = rand(1)*(xMax-xMin) + xMin;
                proposal(2) = rand(1)*(yMax-yMin) + yMin;
                if inShape(aShape,proposal)
                    % Set the proposed point as leaf start point
                    leafSP = proposal;
                    % Increment sampling counter and calculate new ratio
                    rh = proposal(3)/maxHeight;
                    iPCS = find(rh<PCSampling.binEdges,1,'first') - 1;
                    PCSampling.nTotalSampled(iPCS) = ...
                                       PCSampling.nTotalSampled(iPCS) + 1;
                    PCSampling.ratio(iPCS) = PCSampling.nPCSampled(iPCS)...
                                           /PCSampling.nTotalSampled(iPCS);
                    accepted = 1;
                end
            end
        end
    end
    


% Use heightwise and directionwise LADD (or point cloud) in positioning
elseif isempty(fDist_d)

    uniformSampling = false;
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
        iPCS = find(hProposal<PCSampling.binEdges,1,'first') - 1;
        if all([PCSampling.ratio(iPCS) <= PCSampling.weights(iPCS), ...
               uniformSampling == false, ...
               PCSampling.flag == true])
            % Corresponding height index in voxelization
            iz = find(PCSampling.zEdges>hLeaf,1,'first') - 1;
            zSlice = PCSampling.voxels(:,:,iz);
            zSliceCol = zSlice(:);
            zSliceInds = PCSampling.voxelIndex(:,:,iz);
            zSliceIndsCol = zSliceInds(:);
            zSliceCenDir = PCSampling.voxelCenterDirection(:,:,iz);
            zSliceCenDirCol = zSliceCenDir(:);
            % Corresponding voxels around the compass direction
            sectorWidth = 2*pi/10;
            inSector = all([(zSliceCenDirCol < cLeaf+0.5*sectorWidth) ...
                (zSliceCenDirCol > cLeaf-0.5*sectorWidth)], ...
                2);
            sectorSliceIndsCol = zSliceIndsCol(inSector);
            sectorSliceCol = zSliceCol(inSector);
            % If there are no voxels, do uniform sampling
            if sum(sectorSliceCol) == 0
                uniformSampling = true;
                continue
            end
            cumulativeSum = cumsum(sectorSliceCol);
            cumulativePDF = cumulativeSum./cumulativeSum(end);
            k = find(cumulativePDF>rand(1),1,'first');
            xyInds = sectorSliceIndsCol{k};
            ix = xyInds(1);
            iy = xyInds(2);
            % Uniformly sample a point inside the voxel
            proposal = rand(1,3) ...
                       .*[PCSampling.xEdges(ix+1)-PCSampling.xEdges(ix) ...
                          PCSampling.yEdges(iy+1)-PCSampling.yEdges(iy) ...
                          0] ...
                       + [PCSampling.xEdges(ix) PCSampling.yEdges(iy) ...
                          hLeaf];
            leafSP = proposal;
            % Increment sampling counters and calculate new ratio
            PCSampling.nPCSampled(iPCS) = PCSampling.nPCSampled(iPCS) + 1;
            PCSampling.nTotalSampled(iPCS) = ...
                                       PCSampling.nTotalSampled(iPCS) + 1;
            PCSampling.ratio(iPCS) = PCSampling.nPCSampled(iPCS) ...
                                     /PCSampling.nTotalSampled(iPCS);
            accepted = 1;
        end

        if uniformSampling == true
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
            leafSP = (rotation_matrix([0 0 1],cLeaf) ...
                              *[0 1 0]')'*dLeaf ...
                             +[stemCen(1:2) hLeaf];
            % Increment sampling counter and calculate new ratio
            iPCS = find(hLeaf<PCSampling.binEdges,1,'first') - 1;
            PCSampling.nTotalSampled(iPCS) = ...
                                       PCSampling.nTotalSampled(iPCS) + 1;
            PCSampling.ratio(iPCS) = PCSampling.nPCSampled(iPCS)...
                                     /PCSampling.nTotalSampled(iPCS);
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
            leafSP = (rotation_matrix([0 0 1],cLeaf) ...
                              *[0 1 0]')'*dLeaf ...
                             +[stemCen(1:2) hLeaf];
            accepted = 1;
        end
    end
    
end
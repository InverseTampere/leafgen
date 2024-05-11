function [NewCylinderParameters,originalIndex] = preprocess_cylinders(...
                                                    CylinderParameters, ...
                                                    Nodes ...
                                                    )

% Leaf cylinder library minimum and maximum cylinder lenghts
lMin = Nodes.cylinderLength(1);
lMax = Nodes.cylinderLength(2);

% Read cylinder parameters (unsorted)
relativeHeight_US   = CylinderParameters.relative_height;
relativePosition_US = CylinderParameters.relative_position;
branchOrder_US      = CylinderParameters.branch_order;
cylAxis_US          = CylinderParameters.axis;
cylLength_US        = CylinderParameters.length;
cylRadius_US        = CylinderParameters.radius;
isLast_US           = CylinderParameters.is_last;
startPoint_US       = CylinderParameters.start_point;
branchIndex_US      = CylinderParameters.branch_index;
indexInBranch_US    = CylinderParameters.index_in_branch;

% Make sure the cylinders are sorted with respect to branch indices and 
% their relative positions within the branch
relativeHeight   = zeros(size(relativeHeight_US));
relativePosition = zeros(size(relativePosition_US));
branchOrder      = zeros(size(branchOrder_US));
cylAxis          = zeros(size(cylAxis_US));
cylLength        = zeros(size(cylLength_US));
cylRadius        = zeros(size(cylRadius_US));
isLast           = false(size(isLast_US));
startPoint       = zeros(size(startPoint_US));
branchIndex      = zeros(size(branchIndex_US));
indexInBranch    = zeros(size(indexInBranch_US));

ogIndices = zeros(length(relativeHeight_US),1);

k = 1;
for iBranch = 1:max(branchIndex_US)
    % Vector indices of cylinders in branch with current branch index
    indBranch = (branchIndex_US == iBranch);
    % Skip branch if there are no cylinders
    if not(any(indBranch))
        continue
    end
    % The permutation of indices to sorting the cylinders of the branch
    [~,permCyl] = sort(indexInBranch_US(indBranch));

    % Indices of the cylinders in the branch
    indCyl = find(indBranch);

    % Add sorted cylinders to cylinder parameter variables
    kk = k + length(indCyl) - 1;

    relativeHeight(k:kk)   = relativeHeight_US(indCyl(permCyl));
    relativePosition(k:kk) = relativePosition_US(indCyl(permCyl));
    branchOrder(k:kk)      = branchOrder_US(indCyl(permCyl));
    cylAxis(k:kk,:)        = cylAxis_US(indCyl(permCyl),:);
    cylLength(k:kk)        = cylLength_US(indCyl(permCyl));
    cylRadius(k:kk)        = cylRadius_US(indCyl(permCyl));
    isLast(k:kk)           = isLast_US(indCyl(permCyl));
    startPoint(k:kk,:)     = startPoint_US(indCyl(permCyl),:);
    branchIndex(k:kk)      = branchIndex_US(indCyl(permCyl));
    indexInBranch(k:kk)    = indexInBranch_US(indCyl(permCyl));

    % Keep track of the sorting of original indices
    ogIndices(k:kk) = indCyl(permCyl);

    k = kk + 1;
end

% Cylinder midpoints
midPoint = startPoint + 0.5.*cylLength.*cylAxis;

% Maximum and minimum height of a cylinder
maxHeight = max(midPoint(:,3));
minHeight = min(midPoint(:,3));

% Total number of cylinders before preprocessing
nCyl = size(cylLength,1);

% Preallocate new cylinder parameter values
newRelativeHeight   = zeros(2*nCyl,1);
newRelativePosition = zeros(2*nCyl,1);
newBranchOrder      = zeros(2*nCyl,1);
newCylAxis          = zeros(2*nCyl,3);
newCylLength        = zeros(2*nCyl,1);
newCylRadius        = zeros(2*nCyl,1);
newIsLast           = false(2*nCyl,1);
newStartPoint       = zeros(2*nCyl,3);
newBranchIndex      = zeros(2*nCyl,1);
newIndexInBranch    = zeros(2*nCyl,1);

% Counter for additional index in branch caused by the splitting of long
% cylinders or the combination of small cylinders
iibCounter = zeros(max(branchIndex),1);

% Variable for marking the original indices of the cylinders
originalIndex = zeros(2*nCyl,1);

% Loop over all cylinders
newInd = 0;
nExtend = 0;
nSkip = 0;
for iCyl = 1:nCyl

    % Skip cylinders if necessary
    if nSkip > 0
        nSkip = nSkip - 1;
        continue
    end

    % Cut too long cylinders into multiple shorter ones
    if cylLength(iCyl) > lMax
        k = ceil(cylLength(iCyl)/lMax);
        newLen = cylLength(iCyl)/k;
        % Make sure the cylinders are longer than the minimum required
        % lenght, otherwise accept the original cylinder
        if newLen > lMin
            for ii = 1:k
                newInd = newInd + 1;
                % New start point
                sp = startPoint(iCyl,:) + (ii-1)*newLen*cylAxis(iCyl,:);
                % New mid point
                mp = sp + 0.5*newLen*cylAxis(iCyl,:);
                % Mark last cylinder if necessary
                if isLast(iCyl) == true && ii == k
                    il = true;
                else
                    il = false;
                end
                % Find the new relative heigth of cylinder
                rh = (mp(3)-minHeight)/(maxHeight-minHeight);
                % Find the relative position in branch
                if iCyl == 1
                    rp = (ii/k)*relativePosition(iCyl);
                elseif branchIndex(iCyl) ~= branchIndex(iCyl-1)
                    rp = (ii/k)*relativePosition(iCyl);
                else 
                    df = relativePosition(iCyl) - relativePosition(iCyl-1);
                    rp = (ii/k)*df;
                end
                % Set the new cylinder parameter values
                if rh > 1 || rh < 0
                    error("negatiivinen arvo suhteelliselle korkeudelle")
                end
                newRelativeHeight(newInd)   = rh;
                newRelativePosition(newInd) = rp;
                newBranchOrder(newInd)      = branchOrder(iCyl);
                newCylAxis(newInd,:)        = cylAxis(iCyl,:);
                newCylLength(newInd)        = newLen;
                newCylRadius(newInd)        = cylRadius(iCyl);
                newIsLast(newInd)           = il;
                newStartPoint(newInd,:)     = sp;
                newBranchIndex(newInd)      = branchIndex(iCyl);
                newIndexInBranch(newInd) = indexInBranch(iCyl) ...
                                           + (ii-1) ...
                                           + iibCounter(branchIndex(iCyl));
                originalIndex(newInd)    = ogIndices(iCyl);
            end
            % Increment index in branch counter
            bi = branchIndex(iCyl);
            iibCounter(bi) = iibCounter(bi) + (k-1); 
        else
            % Include the cylinder in original length
            newInd = newInd + 1;
            newRelativeHeight(newInd)   = relativeHeight(iCyl);
            newRelativePosition(newInd) = relativePosition(iCyl);
            newBranchOrder(newInd)      = branchOrder(iCyl);
            newCylAxis(newInd,:)        = cylAxis(iCyl,:);
            newCylLength(newInd)        = cylLength(iCyl);
            newCylRadius(newInd)        = cylRadius(iCyl);
            newIsLast(newInd)           = isLast(iCyl);
            newStartPoint(newInd,:)     = startPoint(iCyl,:);
            newBranchIndex(newInd)      = branchIndex(iCyl);
            newIndexInBranch(newInd)    = indexInBranch(iCyl) ...
                                          + iibCounter(branchIndex(iCyl));
            originalIndex(newInd)       = ogIndices(iCyl);
        end


    elseif cylLength(iCyl) < lMin
        % Branch index of the current cylinder
        biCyl = branchIndex(iCyl);
        % Combined length of the included cylinders
        combLen = cylLength(iCyl);
        % Branch index of the next cylinder
        if iCyl < nCyl
            biNext = branchIndex(iCyl+1);
        else
            biNext = -1;
        end
        % Cylinder index counter
        ii = iCyl;
        % Combine cylinders until minimum length is achieved or the
        % cylinders in branch run out
        while combLen < lMin && biNext == biCyl
            ii = ii + 1;
            combLen = combLen + cylLength(ii);
            if ii+1 <= nCyl
                biNext = branchIndex(ii+1);
            else
                break
            end
        end
        newStart = startPoint(iCyl,:);
        newEnd   = startPoint(ii,:) + cylLength(ii)*cylAxis(ii,:);
        newAxis  = (newEnd-newStart)/norm(newEnd-newStart);
        newLength = combLen;
        newRadius = mean(cylRadius(iCyl:ii));
        newMidPoint = newStart + 0.5*newLength*newAxis;
        
        % Set the new cylinder parameter values
        newInd = newInd + 1;
        newRelativeHeight(newInd)   = (newMidPoint(3)-minHeight) ...
                                      /(maxHeight-minHeight);
        newRelativePosition(newInd) = relativePosition(ii);
        newBranchOrder(newInd)      = biCyl;
        newCylAxis(newInd,:)        = newAxis;
        newCylLength(newInd)        = newLength;
        newCylRadius(newInd)        = newRadius;
        newIsLast(newInd)           = isLast(ii);
        newStartPoint(newInd,:)     = startPoint(iCyl,:);
        newBranchIndex(newInd)      = branchIndex(iCyl);
        newIndexInBranch(newInd)    = indexInBranch(iCyl) ...
                                      + iibCounter(branchIndex(iCyl));
        originalIndex(newInd)       = ogIndices(ii);

        % Decrease the index in branch counter
        iibCounter(biCyl) = iCyl - ii;

        % Skip the subsequent combined cylinders
        nSkip = ii - iCyl;

    else
        % Include the cylinder in original length
        newInd = newInd + 1;
        newRelativeHeight(newInd)   = relativeHeight(iCyl);
        newRelativePosition(newInd) = relativePosition(iCyl);
        newBranchOrder(newInd)      = branchOrder(iCyl);
        newCylAxis(newInd,:)        = cylAxis(iCyl,:);
        newCylLength(newInd)        = cylLength(iCyl);
        newCylRadius(newInd)        = cylRadius(iCyl);
        newIsLast(newInd)           = isLast(iCyl);
        newStartPoint(newInd,:)     = startPoint(iCyl,:);
        newBranchIndex(newInd)      = branchIndex(iCyl);
        newIndexInBranch(newInd)    = indexInBranch(iCyl) ...
                                      + iibCounter(branchIndex(iCyl));
        originalIndex(newInd)       = ogIndices(iCyl);

    end

    % If necessary, increase the size of preallocated variables
    if newInd > (nExtend+1)*2*nCyl

        nExtend = nExtend + 1;

        newRelativeHeight   = [newRelativeHeight;   zeros(2*nCyl,1)]; %#ok<AGROW>
        newRelativePosition = [newRelativePosition; zeros(2*nCyl,1)]; %#ok<AGROW>
        newBranchOrder      = [newBranchOrder;      zeros(2*nCyl,1)]; %#ok<AGROW>
        newCylAxis          = [newCylAxis;          zeros(2*nCyl,3)]; %#ok<AGROW>
        newCylLength        = [newCylLength;        zeros(2*nCyl,1)]; %#ok<AGROW>
        newCylRadius        = [newCylRadius;        zeros(2*nCyl,1)]; %#ok<AGROW>
        newIsLast           = [newIsLast;           false(2*nCyl,1)]; %#ok<AGROW>
        newStartPoint       = [newStartPoint;       zeros(2*nCyl,3)]; %#ok<AGROW>
        newBranchIndex      = [newBranchIndex;      zeros(2*nCyl,1)]; %#ok<AGROW>
        newIndexInBranch    = [newIndexInBranch;    zeros(2*nCyl,1)]; %#ok<AGROW>

        originalIndex       = [originalIndex;       zeros(2*nCyl,1)]; %#ok<AGROW>
    end
    
end

originalIndex = originalIndex(1:newInd);

% Return new cylinder parameter values
NewCylinderParameters.relative_height   = newRelativeHeight(1:newInd);
NewCylinderParameters.relative_position = newRelativePosition(1:newInd);
NewCylinderParameters.branch_order      = newBranchOrder(1:newInd);
NewCylinderParameters.axis              = newCylAxis(1:newInd,:);
NewCylinderParameters.length            = newCylLength(1:newInd);
NewCylinderParameters.radius            = newCylRadius(1:newInd);
NewCylinderParameters.is_last           = newIsLast(1:newInd);
NewCylinderParameters.start_point       = newStartPoint(1:newInd,:);
NewCylinderParameters.branch_index      = newBranchIndex(1:newInd);
NewCylinderParameters.index_in_branch   = newIndexInBranch(1:newInd);

end
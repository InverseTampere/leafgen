function RDAS = calculate_relative_distance_along_subbranch( ...
                    CylinderParameters)
% Set the reference point of cylinder position
referencePoint = "end";

% Branch indexes
branchIndex = CylinderParameters.branch_index;
% Lenghts of cylinders
cylinderLength = CylinderParameters.length;

% Total number of cylinders
nCyl = length(branchIndex);

% Relative along-branch distances from base of sub-branch
RDAS = zeros(nCyl,1);
indexVector = (1:1:nCyl)';
for iBranch = 1:max(branchIndex) % stem index 0 is skipped automatically
    % Indexes of QSM cylinders belonging in branch
    bcIndexes = indexVector(branchIndex == iBranch);
    if sum(bcIndexes) == 0
        % the branch has no cylinders
        continue
    end
    bcIndexes = nonzeros(bcIndexes);
    % Lengths of corresponding cylinders
    bcLengths = zeros(length(bcIndexes),1);
    for j = 1:length(bcLengths)
        bcLengths(j) = cylinderLength(bcIndexes(j));
    end
    % Cumulative sum of the lengths
    bcLengthsCumulative = cumsum(bcLengths);
    if referencePoint == "end"
        % Relative position in subbranch
        relDistInSubBranch = bcLengthsCumulative/bcLengthsCumulative(end);
    elseif referencePoint == "middle"
        % Cumulative distance of cylinder midpoints in sub-branch
        midPointCumulative = bcLengthsCumulative - 0.5*bcLengths;
        % Relative position in subbranch
        relDistInSubBranch = midPointCumulative/bcLengthsCumulative(end);
    end
    k = 1;
    for j = 1:length(bcIndexes)
        bcInd = bcIndexes(j);
        RDAS(bcInd) = relDistInSubBranch(k);
        k = k + 1;
    end
end
end
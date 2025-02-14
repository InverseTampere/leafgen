function [leafDir,leafNormal,petioleStart,petioleEnd] = fun_leaf_orientation( ...
                                                    len,rad,inc,az, ...
                                                    nLeaves, ...
                                                    petioleLengthLimits,...
                                                    dTypeLODinc, ...
                                                    dTypeLODaz, ...
                                                    dParametersLODinc, ...
                                                    dParametersLODaz, ...
                                                    PetioleDirDistribution,...
                                                    Phyllotaxis)

% Function definitions
fun_beta = @(x,a,b) (1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
fun_vonmises = @(x,m,k) exp(k*cos(x-m))./(2*pi*besseli(0,k));
fun_dewit = @(x,a,b) (1 + a*cos(b*x))/(pi/2+(a/b)*sin(b*pi/2));

% Define the distribution function for leaf inclination angle
switch dTypeLODinc

    case 'uniform'
        % Uniform distribution
        f_inc = @(x,p) 2/pi*ones(size(x));
        max_f_inc = @(p) 2/pi;
        % Set sampling type to rejection sampling
        lod_inc_sampling = 'rejection sampling';

    case 'spherical'
        % Spherical distribution function
        f_inc = @(x,p) sin(x);
        % Inverse of cumulative density function
        F_inc_inv = @(y,p) acos(1-y);
        % Set sampling type to inverse sampling
        lod_inc_sampling = 'inverse sampling';

    case 'dewit'
        % Generalized de Wit's distribution function
        f_inc = @(x,p) fun_dewit(x,p(1),p(2));
        % Upper limit for the distribution
        max_f_inc = @(p) max(f_inc([0 1 2 3]*(pi/p(2)),p));
        % Set sampling type to rejection sampling
        lod_inc_sampling = 'rejection sampling';

    case 'beta'
        % Beta distribution density function
        f_inc = @(x,p) fun_beta(2*x/pi,p(1),p(2));
        % Inverse of cumulative density function
        F_inc_inv = @(y,p) (pi/2)*betaincinv(y,p(1),p(2));
        % Set sampling type to inverse sampling
        lod_inc_sampling = 'inverse sampling';

    case 'constant'
        % Set sampling type to constant
        lod_inc_sampling = 'constant';

end

% Define the distribution function for leaf azimuth angle
switch dTypeLODaz

    case 'uniform'
        % Uniform distribution
        f_az = @(x,p) 1/(2*pi)*ones(size(x));
        max_f_az = @(p) 1/(2*pi);
        % Set sampling type to rejection sampling
        lod_az_sampling = 'rejection sampling';

    case 'vonmises'
        % Von Mises distribution density function
        f_az = @(x,p) fun_vonmises(x,p(1),p(2));
        % Upper limit for the distribution
        max_f_az = @(p) fun_vonmises(p(1),p(1),p(2));
        % Set sampling type to rejection sampling
        lod_az_sampling = 'rejection sampling';

    case 'constant'
        % Set sampling type to constant
        lod_az_sampling = 'constant';

end

% Unit axis vector of the cylinder
cylinderAxis = [sin(inc)*cos(az+pi/2), ...
                sin(inc)*sin(az+pi/2), ...
                cos(inc)];

incAngles = zeros(nLeaves,1);
azAngles  = zeros(nLeaves,1);
leafNormal = zeros(nLeaves,3);
leafDir = zeros(nLeaves,3);
petioleStart = zeros(nLeaves,3);
petioleEnd = zeros(nLeaves,3);

if Phyllotaxis.flag == true
    % Relative lengthwise node position on cylinder axis (starts from the
    % end tip of the cylinder)
    nodeRelPosAxis = 1;
    % Set the radial direction angle value for the first petiole
    if isfield(Phyllotaxis,'initialRadialAngle')
        nodeRotationAngle = Phyllotaxis.initialRadialAngle;
    else
        nodeRotationAngle = 0;
    end
end

% Variable for skipping leaf indices in case of phyllotaxis
skipIndices = 0;

for iLeaf = 1:nLeaves

    % Skip indices if phyllotaxis pattern has multiple leaves per node
    if skipIndices > 0
        skipIndices = skipIndices - 1;
        continue
    end
        

    % Sample the leaf inclination angle
    switch lod_inc_sampling

        case 'rejection sampling'
            % Sample inclination angle value with acceptance-rejection
            % sampling
            accepted = 0;
            while accepted == 0
                incProposal = rand(1)*pi/2;
                funValue = f_inc(incProposal,dParametersLODinc);
                vertValue = rand(1)*max_f_inc(dParametersLODinc);
                if vertValue < funValue
                    incAngles(iLeaf) = incProposal;
                    accepted = 1;
                end
            end

        case 'inverse sampling'
            % Sample inclination angle by inverse transform sampling
            u = rand(1);
            incAngles(iLeaf) = F_inc_inv(u,dParametersLODinc);

        case 'constant'
            incAngles(iLeaf) = dParametersLODinc(1);

    end

    % Sample the leaf azimuth angle
    switch lod_az_sampling

        case 'rejection sampling'
            % Sample azimuth angle value with acceptance-rejection
            % sampling
            accepted = 0;
            while accepted == 0
                azProposal = rand(1)*2*pi;
                funValue = f_az(azProposal,dParametersLODaz);
                vertValue = rand(1)*max_f_az(dParametersLODaz);
                if vertValue < funValue
                    azAngles(iLeaf) = azProposal;
                    accepted = 1;
                end
            end

        case 'constant'
            azAngles(iLeaf) = dParametersLODaz(1);
            
    end

    % Unit vector of the leaf normal (y-axis assumed as north direction)
    leafNormal(iLeaf,:) = [sin(incAngles(iLeaf))*cos(azAngles(iLeaf)+pi/2),...
                           sin(incAngles(iLeaf))*sin(azAngles(iLeaf)+pi/2),...
                           cos(incAngles(iLeaf))];

    % Define petiole and leaf direction
    if isfield(Phyllotaxis,'nodeRotation')
        phi = nodeRotationAngle;
    elseif PetioleDirDistribution.flag == true % user-supplied distribution
        f_petiole_dir = PetioleDirDistribution.dist_fun;
        maxVal = max(f_petiole_dir(0:0.001:2*pi,len,rad,inc,az)) + 0.1;
        accepted = 0;
        while accepted == 0
            proposal = rand(1)*2*pi;
            funValue = f_petiole_dir(proposal,len,rad,inc,az);
            vertValue = rand(1)*maxVal;
            if vertValue < funValue
                phi = proposal;
                accepted = 1;
            end
        end
    else
        phi = 2*pi*rand(1);
    end
    
    % If inclination angle is nonzero, choose the upmost direction on
    % the plane defined by leaf normal as the reference direction
    if inc > 1e-3 && inc < (pi - 1e-3)
        % Cylinder side direction unit vector
        cylSide = cross(cylinderAxis,[0 0 1]);
        cylSide = cylSide/norm(cylSide);
        % Upmost pointing radial direction unit vector
        upmostRadDir = cross(cylSide,cylinderAxis);
        upmostRadDir = upmostRadDir/norm(upmostRadDir);
        if upmostRadDir(3) < 0
            upmostRadDir = -upmostRadDir;
        end
        % Find petiole start radial direction by rotating upmost direction 
        % wrt. cylinder axis
        petioleStartRadDir = rotation_matrix(cylinderAxis,phi)*upmostRadDir';
        petioleStartRadDir = petioleStartRadDir';

    % If inclination angle is close to zero, choose the northmost
    % direction on the plane defined by leaf normal as the reference
    % direction
    else
        if cylinderAxis(3) > 0
            northmostDir = cross([-1 0 0],cylinderAxis);
            northmostDir = northmostDir/norm(northmostDir);
            petioleStartRadDir = rotation_matrix(cylinderAxis,phi) ...
                              *northmostDir';
        else
            northmostDir = cross([ 1 0 0],cylinderAxis);
            northmostDir = northmostDir/norm(northmostDir);
            petioleStartRadDir = rotation_matrix(-cylinderAxis,phi) ...
                              *northmostDir';
        end
        petioleStartRadDir = petioleStartRadDir';
    end
    
    % Side direction of petiole start point
    if norm(petioleStartRadDir-leafNormal(iLeaf,:)) > 1e-3 && ...
       norm(petioleStartRadDir+leafNormal(iLeaf,:)) > 1e-3 % cross prod ok
        petioleStartSide = cross(petioleStartRadDir,leafNormal(iLeaf,:));
    else % radial direction of petiole start parallel to leaf normal
        % Make a small perturbation to radial petiole start direction
        petioleStartRadDir = rotation_matrix(cylinderAxis,pi/6-(rand(1)*pi/3)) ...
                          *petioleStartRadDir';
        petioleStartRadDir = petioleStartRadDir';
        petioleStartSide = cross(petioleStartRadDir,leafNormal(iLeaf,:));
    end
    petioleStartSide = petioleStartSide/norm(petioleStartSide);
    
    % Petiole and leaf direction
    if Phyllotaxis.flag == true && isfield(Phyllotaxis, ...
                                           'petioleAxialInclinationAngle')
        rotAxis = cross(cylinderAxis,petioleStartRadDir);
        rotAxis = rotAxis/norm(rotAxis);
        rmPetiole = rotation_matrix(rotAxis, ...
                                 Phyllotaxis.petioleAxialInclinationAngle);
        petioleDir =  (rmPetiole*cylinderAxis')';
        leafDir(iLeaf,:) = cross(leafNormal(iLeaf,:),petioleStartSide);
        leafDir = leafDir/norm(leafDir);
    else
        petioleDir = cross(leafNormal(iLeaf,:),petioleStartSide);
        petioleDir = petioleDir/norm(petioleDir);
        % Set leaf direction to be the same as petiole direction
        leafDir(iLeaf,:) = petioleDir;
    end
    
    % Lengthwise petiole start location on cylinder axis
    if Phyllotaxis.flag == true
        petioleStartPosAxis = nodeRelPosAxis*len;
        nodeRelPosAxis = (nodeRelPosAxis*len-Phyllotaxis.nodeDistance)/len;
    else
        petioleStartPosAxis = rand(1)*len;
    end
    % Start and end points of the petiole
    petioleStart(iLeaf,:) = petioleStartPosAxis*cylinderAxis ...
                         + rad*petioleStartRadDir;
    petioleLen = (petioleLengthLimits(2)-petioleLengthLimits(1))*rand(1) ... 
              + petioleLengthLimits(1);
    petioleEnd(iLeaf,:) = petioleStart(iLeaf,:) + petioleLen*petioleDir;

    % If the phyllotaxis pattern has multiple leaves per node, add them
    if Phyllotaxis.flag == true && Phyllotaxis.nNodeLeaves > 1
        for iPhyl = 1:Phyllotaxis.nNodeLeaves-1
            leafNormal(iLeaf+iPhyl,:) = leafNormal(iLeaf,:);
            rmPhyl = rotation_matrix(cylinderAxis, ...
                                     Phyllotaxis.petioleSeparationAngle);
            petioleStartRadDir = (rmPhyl*petioleStartRadDir')';
            % Side direction of petiole start point
            if norm(petioleStartRadDir-leafNormal(iLeaf,:)) > 1e-3 && ...
                    norm(petioleStartRadDir+leafNormal(iLeaf,:)) > 1e-3 
                petioleStartSide = cross(petioleStartRadDir,leafNormal(iLeaf,:));
            else
                % Make a small perturbation to radial petiole start 
                % direction
                prtb = pi/6-(rand(1)*pi/3);
                petioleStartRadDir = (rotation_matrix(cylinderAxis,prtb) ...
                                   *petioleStartRadDir')';
                petioleStartSide = cross(petioleStartRadDir,leafNormal(iLeaf,:));
            end
            petioleStartSide = petioleStartSide/norm(petioleStartSide);
            % Petiole direction
            if Phyllotaxis.flag == true && ...
               isfield(Phyllotaxis,'petioleAxialInclinationAngle')
                rotAxis = cross(cylinderAxis,petioleStartRadDir);
                rotAxis = rotAxis/norm(rotAxis);
                rmPetiole = rotation_matrix(rotAxis, ...
                                Phyllotaxis.petioleAxialInclinationAngle);
                petioleDir =  (rmPetiole*cylinderAxis')';
            else
                petioleDir = cross(leafNormal(iLeaf,:),petioleStartSide);
                petioleDir = petioleDir/norm(petioleDir);
            end
            % Set leaf direction to be the same as petiole direction
            leafDir(iLeaf+iPhyl,:) = petioleDir;
            % Start and end points of the petiole
            petioleStart(iLeaf+iPhyl,:) = petioleStartPosAxis*cylinderAxis ...
                                       + rad*petioleStartRadDir;
            petioleLen = (petioleLengthLimits(2)-petioleLengthLimits(1))*rand(1) ...
                      + petioleLengthLimits(1);
            petioleEnd(iLeaf+iPhyl,:) = petioleStart(iLeaf+iPhyl,:) ...
                                     + petioleLen*petioleDir;
        end
        skipIndices = Phyllotaxis.nNodeLeaves-1;
    end

    if Phyllotaxis.flag == true
        % Break loop if there is no more room for phyllotaxis pattern
        if nodeRelPosAxis < 0
            break
        end
        % Set phyllotaxis rotation angle of successive nodes
        if isfield(Phyllotaxis,'nodeRotation')
            nodeRotationAngle = nodeRotationAngle ...
                                + Phyllotaxis.nodeRotation;
        end
    end
end

if Phyllotaxis.flag == true
    % Remove empty rows if necessary
    totalAddedLeaves = iLeaf + Phyllotaxis.nNodeLeaves - 1;
    if totalAddedLeaves < nLeaves
        leafNormal = leafNormal(1:totalAddedLeaves,:);
        leafDir = leafDir(1:totalAddedLeaves,:);
        petioleStart = petioleStart(1:totalAddedLeaves,:);
        petioleEnd = petioleEnd(1:totalAddedLeaves,:);
    end
end

end
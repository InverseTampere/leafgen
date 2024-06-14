function QSM = check_inputs_direct(QSM,TargetDistributions, ...
                                   LeafProperties,totalLeafArea);

%% QSM
% Convert QSM to QSMBCylindrical-class object if necessary
if class(QSM) ~= "QSMBCylindrical"
    QSM = QSMBCylindrical(QSM);
end

%% Target distributions
% Assert the existence of all fields of the struct
fieldCheckTar = [isfield(TargetDistributions,'dTypeLADDh'), ...
                 isfield(TargetDistributions,'dTypeLADDd'), ...
                 isfield(TargetDistributions,'dTypeLADDc'), ...
                 isfield(TargetDistributions,'dTypeLODinc'), ...
                 isfield(TargetDistributions,'dTypeLODaz'), ...
                 isfield(TargetDistributions,'dTypeLSD'), ...
                 isfield(TargetDistributions,'hParams'), ...
                 isfield(TargetDistributions,'dParams'), ...
                 isfield(TargetDistributions,'cParams'), ...
                 isfield(TargetDistributions,'fun_inc_params'), ...
                 isfield(TargetDistributions,'fun_az_params'), ...
                 isfield(TargetDistributions,'fun_size_params'), ...
                 ];
fieldNamesTar = ["dTypeLADDh", ...
                 "dTypeLADDd", ...
                 "dTypeLADDc", ...
                 "dTypeLODinc", ...
                 "dTypeLODaz", ...
                 "dTypeLSD", ...
                 "hParams", ...
                 "dParams", ...
                 "cParams", ...
                 "fun_inc_params", ...
                 "fun_az_params", ...
                 "fun_size_params", ...
                  ];
for iField = 1:length(fieldCheckTar)
    assert(fieldCheckTar(iField),"TargetDistributions."...
           +fieldNamesTar(iField)+" is missing.")
end
% Relative height
dTypeH  = TargetDistributions.dTypeLADDh;
hParams = TargetDistributions.hParams;
if ~any(strcmp(dTypeH,{'uniform','polynomial', ...
        'polynomialmixturemodel','weibull','weibullmixturemodel', ...
        'beta','betamixturemodel'}))
    error("LADD height distribution type not recognized.")
end
switch dTypeH
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'polynomial'
        % Assure that the polynomial gets only nonnegative values
        if any(polyval(hParams,0:0.001:1) < 0)
            error("TargetDistributions relative height polynomial gets"...
                  +" negative values on the interval [0,1].")
        end
    case 'polynomialmixturemodel'
        % Pick polynomial coefficients and weight
        nP = (length(hParams)-1)/2; % number of polynomial coefficients
        p1 = hParams(1:nP); % coefficients of the first polynomial
        p2 = hParams((nP+1):(2*nP)); % coefficients of the second polynom.
        w = hParams(end); % mixture model weight
        xx = 0:0.001:1;
        mmPolyValues = w*polyval(p1,xx) ...
                       + (1-w)*polyval(p2,xx);
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.hParams mixture model weight is"...
                  +" not on the interval [0,1].")
        end
        % Assure that the polynomial gets only nonnegative values
        if any(mmPolyValues < 0)
            error("TargetDistributions relative height polynomial gets"...
                  +" negative values on the interval [0,1].")
        end
    case 'weibull'
        % Assure that both parameters are positive
        if any(hParams <= 0)
            error("TargetDistributions.hParams all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
    case 'weibullmixturemodel'
        % Assure that all Weibull parameters are positive
        if any(hParams(1:4) <= 0)
            error("TargetDistributions.hParams all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
        w = hParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.hParams mixture model weight is"...
                  +" not on the interval [0,1].")
        end
    case 'beta'
        % Assure that both parameters are positive
        if any(hParams <= 0)
            error("TargetDistributions.hParams all elements have to be"...
                  +" positive for beta distribution.")
        end
    case 'betamixturemodel'
        % Assure that all beta parameters are positive
        if any(hParams(1:4) <= 0)
            error("TargetDistributions.hParams all elements have to be"...
                  +" positive for beta distribution.")
        end
        w = hParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.hParams mixture model weight is"...
                  +" not on the interval [0,1].")
        end
end
% Relative distance along subbranch
dTypeD  = TargetDistributions.dTypeLADDd;
dParams = TargetDistributions.dParams;
if ~any(strcmp(dTypeD,{'uniform','polynomial', ...
        'polynomialmixturemodel','weibull','weibullmixturemodel', ...
        'beta','betamixturemodel'}))
    error("LADD distance from stem distribution type not recognized.")
end
switch dTypeD
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'polynomial'
        % Assure that the polynomial gets only nonnegative values
        if any(polyval(dParams,0:0.001:1) < 0)
            error("TargetDistributions relative distance along"...
                  +" subbranch polynomial gets negative values on the"...
                  +" interval [0,1].")
        end
    case 'polynomialmixturemodel'
        % Pick polynomial coefficients and weight
        nP = (length(dParams)-1)/2; % number of polynomial coefficients
        p1 = dParams(1:nP); % coefficients of the first polynomial
        p2 = dParams((nP+1):(2*nP)); % coefficients of the second polynom.
        w = dParams(end); % mixture model weight
        xx = 0:0.001:1;
        mmPolyValues = w*polyval(p1,xx) ...
                       + (1-w)*polyval(p2,xx);
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.dParams mixture model weight is"...
                  +" not on the interval [0,1].")
        end
        % Assure that the polynomial gets only nonnegative values
        if any(mmPolyValues < 0)
            error("TargetDistributions relative distance along"...
                  +" subbranch polynomial gets negative values on the"...
                  +" interval [0,1].")
        end
    case 'weibull'
        % Assure that both parameters are positive
        if any(dParams <= 0)
            error("TargetDistributions.dParams all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
    case 'weibullmixturemodel'
        % Assure that all Weibull parameters are positive
        if any(dParams(1:4) <= 0)
            error("TargetDistributions.dParams all elements have to be"...
                  +" positive for truncated Weibull distributions.")
        end
        w = dParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.dParams mixture model weight is"...
                  +" not on the interval [0,1].")
        end
    case 'beta'
        % Assure that both parameters are positive
        if any(dParams <= 0)
            error("TargetDistributions.dParams all elements have to be"...
                  +" positive for beta distribution.")
        end
    case 'betamixturemodel'
        % Assure that all beta parameters are positive
        if any(dParams(1:4) <= 0)
            error("TargetDistributions.dParams all elements have to be"...
                  +" positive for beta distributions.")
        end
        w = dParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.dParams mixture model weight is"...
                  +" not on the interval [0,1].")
        end
end
% Compass direction
dTypeC  = TargetDistributions.dTypeLADDc;
cParams = TargetDistributions.cParams;
if ~any(strcmp(dTypeC,{'uniform','vonmises','vonmisesmixturemodel'}))
    error("LADD compass direction distribution type not recognized.")
end
switch dTypeC
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'vonmises'
        % Assure that the second parameter is positive
        if cParams(2) <= 0
            error("TargetDistributions.cParams second parameter has"...
                  +" to be positive for von Mises distribution.")
        end
    case 'vonmisesmixturemodel'
        % Assure that second parameters are positive for both disrtibutions
        if cParams(2) <= 0 || cParams(4) <= 0
            error("TargetDistributions.cParams second and fourth"...
                  +" element have to be positive for von Mises"...
                  +" distributions.")
        end
        w = cParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.cParams mixture model weight"...
                  +" is not on the interval [0,1].")
        end
end
% Parameter functions
% Assert that all fields are function handles
if ~isa(TargetDistributions.fun_inc_params,'function_handle')
    error("TargetDistributions.fun_inc_params has to be a function"...
          +" handle.")
end
if ~isa(TargetDistributions.fun_az_params,'function_handle')
    error("TargetDistributions.fun_az_params has to be a function handle.")
end
if ~isa(TargetDistributions.fun_size_params,'function_handle')
    error("TargetDistributions.fun_size_params has to be a function"...
          +" handle.")
end

%% Check the leaf base geometry and twig length limits
if size(LeafProperties.vertices,2) ~= 3
    error("LeafProperties.vertices should have a size (nVertices x 3).")
end
if size(LeafProperties.triangles,2) ~= 3
    error("LeafProperties.triangles should have a size (nTriangles x 3).")
end
if max(LeafProperties.triangles,[],'all') > size(LeafProperties.vertices,1)
    error("LeafProperties.triangles refers to too large vertice indices.")
end
if any(LeafProperties.twigLengthLimits < 0) || ...
        ~all(LeafProperties.twigLengthLimits ...
             == sort(LeafProperties.twigLengthLimits))
    error("LeafProperties.twigLengthLimits can contain only positive"...
          +" values in ascending order.")
end
if length(LeafProperties.twigLengthLimits) ~= 2
    error("LeafProperties.twigLengthLimits can contain only two elements.")
end

%% Total leaf area
% Assert that total leaf area is positive
if totalLeafArea <= 0
    error("totalLeafArea has to be positive.")
end

end
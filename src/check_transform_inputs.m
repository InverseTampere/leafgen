function QSM = check_transform_inputs(QSM,TargetLADD,ParamFunctions, ...
                                      targetLeafArea)

%% QSM
% Convert QSM to QSMBCylindrical-class object if necessary
if class(QSM) ~= "QSMBCylindrical"
    QSM = QSMBCylindrical(QSM);
end

%% Target LADD
% Assert the existence of all fields of the struct
fieldCheckLADD = [isfield(TargetLADD,'dTypeLADDh'), ...
                  isfield(TargetLADD,'dTypeLADDd'), ...
                  isfield(TargetLADD,'dTypeLADDc'), ...
                  isfield(TargetLADD,'hParams'), ...
                  isfield(TargetLADD,'dParams'), ...
                  isfield(TargetLADD,'cParams'), ...
                  ];
fieldNamesLADD = ["dTypeLADDh", ...
                  "dTypeLADDd", ...
                  "dTypeLADDc", ...
                  "hParams", ...
                  "dParams", ...
                  "cParams", ...
                  ];
for iField = 1:length(fieldCheckLADD)
    assert(fieldCheckLADD(iField),"TargetLADD."+fieldNamesLADD(iField)+...
           " is missing.")
end
% Relative height
dTypeH  = TargetLADD.dTypeLADDh;
hParams = TargetLADD.hParams;
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
            error("TargetLADD relative height polynomial gets negative"...
                  +" values on the interval [0,1].")
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
            error("TargetLADD.hParams mixture model weight is not on"...
                  +" the interval [0,1].")
        end
        % Assure that the polynomial gets only nonnegative values
        if any(mmPolyValues < 0)
            error("TargetLADD relative height polynomial gets negative"...
                  +" values on the interval [0,1].")
        end
    case 'weibull'
        % Assure that both parameters are positive
        if any(hParams <= 0)
            error("TargetLADD.hParams all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
    case 'weibullmixturemodel'
        % Assure that all Weibull parameters are positive
        if any(hParams(1:4) <= 0)
            error("TargetLADD.hParams all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
        w = hParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.hParams mixture model weight is not on"...
                  +" the interval [0,1].")
        end
    case 'beta'
        % Assure that both parameters are positive
        if any(hParams <= 0)
            error("TargetLADD.hParams all elements have to be"...
                  +" positive for beta distribution.")
        end
    case 'betamixturemodel'
        % Assure that all beta parameters are positive
        if any(hParams(1:4) <= 0)
            error("TargetLADD.hParams all elements have to be"...
                  +" positive for beta distribution.")
        end
        w = hParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.hParams mixture model weight is not on"...
                  +" the interval [0,1].")
        end
end
% Relative distance along subbranch
dTypeD = TargetLADD.dTypeLADDd;
dParams = TargetLADD.dParams;
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
            error("TargetLADD relative distance along subbranch"...
                  +" polynomial gets negative values on the interval"...
                  +" [0,1].")
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
            error("TargetLADD.dParams mixture model weight is not on"...
                  +" the interval [0,1].")
        end
        % Assure that the polynomial gets only nonnegative values
        if any(mmPolyValues < 0)
            error("TargetLADD relative distance along subbranch"...
                  +" polynomial gets negative values on the interval"...
                  +" [0,1].")
        end
    case 'weibull'
        % Assure that both parameters are positive
        if any(dParams <= 0)
            error("TargetLADD.dParams all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
    case 'weibullmixturemodel'
        % Assure that all Weibull parameters are positive
        if any(dParams(1:4) <= 0)
            error("TargetLADD.dParams all elements have to be"...
                  +" positive for truncated Weibull distributions.")
        end
        w = dParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.dParams mixture model weight is not on"...
                  +" the interval [0,1].")
        end
    case 'beta'
        % Assure that both parameters are positive
        if any(dParams <= 0)
            error("TargetLADD.dParams all elements have to be"...
                  +" positive for beta distribution.")
        end
    case 'betamixturemodel'
        % Assure that all beta parameters are positive
        if any(dParams(1:4) <= 0)
            error("TargetLADD.dParams all elements have to be"...
                  +" positive for beta distributions.")
        end
        w = dParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.dParams mixture model weight is not on"...
                  +" the interval [0,1].")
        end
end
% Compass direction
dTypeC  = TargetLADD.dTypeLADDc;
cParams = TargetLADD.cParams;
if ~any(strcmp(dTypeC,{'uniform','vonmises','vonmisesmixturemodel'}))
    error("LADD compass direction distribution type not recognized.")
end
switch dTypeC
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'vonmises'
        % Assure that the second parameter is positive
        if cParams(2) <= 0
            error("TargetLADD.cParams second parameter has to be"...
                  +" positive for von Mises distribution.")
        end
    case 'vonmisesmixturemodel'
        % Assure that second parameters are positive for both disrtibutions
        if cParams(2) <= 0 || cParams(4) <= 0
            error("TargetLADD.cParams second and fourth element have"...
                  +" to be positive for von Mises distributions.")
        end
        w = cParams(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.cParams mixture model weight is not on"...
                  +" the interval [0,1].")
        end
end

%% Parameter functions
% Assert the existence of all fields of the struct
fieldCheckPar = [isfield(ParamFunctions,'fun_inc_params'), ...
                 isfield(ParamFunctions,'fun_az_params'), ...
                 isfield(ParamFunctions,'fun_size_params'), ...
                 ];
fieldNamesPar = ["fun_inc_params", ...
                 "fun_az_params", ...
                 "fun_size_params", ...
                 ];
for iField = 1:length(fieldCheckPar)
    assert(fieldCheckPar(iField),"ParamFunctoins."+...
           fieldNamesPar(iField)+" is missing.")
end
% Assert that all field are function handles
if ~isa(ParamFunctions.fun_inc_params,'function_handle')
    error("ParamFunctions.fun_inc_params has to be a function handle.")
end
if ~isa(ParamFunctions.fun_az_params,'function_handle')
    error("ParamFunctions.fun_az_params has to be a function handle.")
end
if ~isa(ParamFunctions.fun_size_params,'function_handle')
    error("ParamFunctions.fun_size_params has to be a function handle.")
end

%% Target leaf area
% Assert that target leaf area is positive
if targetLeafArea <= 0
    error("targetLeafArea has to be positive.")
end
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

function QSM = check_inputs_transform(QSM,TargetLADD,ParamFunctions, ...
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
                  isfield(TargetLADD,'pLADDh'), ...
                  isfield(TargetLADD,'pLADDd'), ...
                  isfield(TargetLADD,'pLADDc'), ...
                  ];
fieldNamesLADD = ["dTypeLADDh", ...
                  "dTypeLADDd", ...
                  "dTypeLADDc", ...
                  "pLADDh", ...
                  "pLADDd", ...
                  "pLADDc", ...
                  ];
for iField = 1:length(fieldCheckLADD)
    assert(fieldCheckLADD(iField),"TargetLADD."+fieldNamesLADD(iField)+...
           " is missing.")
end
% Relative height
dTypeH  = TargetLADD.dTypeLADDh;
pLADDh = TargetLADD.pLADDh;
if ~any(strcmp(dTypeH,{'uniform','polynomial','polynomialmixture', ...
        'weibull','weibullmixture','beta','betamixture'}))
    error("LADD height distribution type not recognized.")
end
switch dTypeH
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'polynomial'
        % Assure that the polynomial gets only nonnegative values
        if any(polyval(pLADDh,0:0.001:1) < 0)
            error("TargetLADD relative height polynomial gets negative"...
                  +" values on the interval [0,1].")
        end
    case 'polynomialmixture'
        % Pick polynomial coefficients and weight
        nP = (length(pLADDh)-1)/2; % number of polynomial coefficients
        p1 = pLADDh(1:nP); % coefficients of the first polynomial
        p2 = pLADDh((nP+1):(2*nP)); % coefficients of the second polynom.
        w = pLADDh(end); % mixture model weight
        xx = 0:0.001:1;
        mmPolyValues = w*polyval(p1,xx) ...
                       + (1-w)*polyval(p2,xx);
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.pLADDh mixture model weight is not on"...
                  +" the interval [0,1].")
        end
        % Assure that the polynomial gets only nonnegative values
        if any(mmPolyValues < 0)
            error("TargetLADD relative height polynomial gets negative"...
                  +" values on the interval [0,1].")
        end
    case 'weibull'
        % Assure that both parameters are positive
        if any(pLADDh <= 0)
            error("TargetLADD.pLADDh all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
    case 'weibullmixture'
        % Assure that all Weibull parameters are positive
        if any(pLADDh(1:4) <= 0)
            error("TargetLADD.pLADDh all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
        w = pLADDh(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.pLADDh mixture model weight is not on"...
                  +" the interval [0,1].")
        end
    case 'beta'
        % Assure that both parameters are positive
        if any(pLADDh <= 0)
            error("TargetLADD.pLADDh all elements have to be"...
                  +" positive for beta distribution.")
        end
    case 'betamixture'
        % Assure that all beta parameters are positive
        if any(pLADDh(1:4) <= 0)
            error("TargetLADD.pLADDh all elements have to be"...
                  +" positive for beta distribution.")
        end
        w = pLADDh(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.pLADDh mixture model weight is not on"...
                  +" the interval [0,1].")
        end
end
% Relative distance along subbranch
dTypeD  = TargetLADD.dTypeLADDd;
pLADDd = TargetLADD.pLADDd;
if ~any(strcmp(dTypeD,{'uniform','polynomial','polynomialmixture', ...
        'weibull','weibullmixture','beta','betamixture'}))
    error("LADD distance from stem distribution type not recognized.")
end
switch dTypeD
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'polynomial'
        % Assure that the polynomial gets only nonnegative values
        if any(polyval(pLADDd,0:0.001:1) < 0)
            error("TargetLADD relative distance along subbranch"...
                  +" polynomial gets negative values on the interval"...
                  +" [0,1].")
        end
    case 'polynomialmixture'
        % Pick polynomial coefficients and weight
        nP = (length(pLADDd)-1)/2; % number of polynomial coefficients
        p1 = pLADDd(1:nP); % coefficients of the first polynomial
        p2 = pLADDd((nP+1):(2*nP)); % coefficients of the second polynom.
        w = pLADDd(end); % mixture model weight
        xx = 0:0.001:1;
        mmPolyValues = w*polyval(p1,xx) ...
                       + (1-w)*polyval(p2,xx);
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.pLADDd mixture model weight is not on"...
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
        if any(pLADDd <= 0)
            error("TargetLADD.pLADDd all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
    case 'weibullmixture'
        % Assure that all Weibull parameters are positive
        if any(pLADDd(1:4) <= 0)
            error("TargetLADD.pLADDd all elements have to be"...
                  +" positive for truncated Weibull distributions.")
        end
        w = pLADDd(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.pLADDd mixture model weight is not on"...
                  +" the interval [0,1].")
        end
    case 'beta'
        % Assure that both parameters are positive
        if any(pLADDd <= 0)
            error("TargetLADD.pLADDd all elements have to be"...
                  +" positive for beta distribution.")
        end
    case 'betamixture'
        % Assure that all beta parameters are positive
        if any(pLADDd(1:4) <= 0)
            error("TargetLADD.pLADDd all elements have to be"...
                  +" positive for beta distributions.")
        end
        w = pLADDd(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.pLADDd mixture model weight is not on"...
                  +" the interval [0,1].")
        end
end
% Compass direction
dTypeC  = TargetLADD.dTypeLADDc;
pLADDc = TargetLADD.pLADDc;
if ~any(strcmp(dTypeC,{'uniform','vonmises','vonmisesmixture'}))
    error("LADD compass direction distribution type not recognized.")
end
switch dTypeC
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'vonmises'
        % Assure that the second parameter is positive
        if pLADDc(2) <= 0
            error("TargetLADD.pLADDc second parameter has to be"...
                  +" positive for von Mises distribution.")
        end
    case 'vonmisesmixture'
        % Assure that second parameters are positive for both disrtibutions
        if pLADDc(2) <= 0 || pLADDc(4) <= 0
            error("TargetLADD.pLADDc second and fourth element have"...
                  +" to be positive for von Mises distributions.")
        end
        w = pLADDc(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetLADD.pLADDc mixture model weight is not on"...
                  +" the interval [0,1].")
        end
end

%% Parameter functions
% Assert the existence of all fields of the struct
fieldCheckPar = [isfield(ParamFunctions,'fun_pLODinc'), ...
                 isfield(ParamFunctions,'fun_pLODaz'), ...
                 isfield(ParamFunctions,'fun_pLSD'), ...
                 ];
fieldNamesPar = ["fun_pLODinc", ...
                 "fun_pLODaz", ...
                 "fun_pLSD", ...
                 ];
for iField = 1:length(fieldCheckPar)
    assert(fieldCheckPar(iField),"ParamFunctions."+...
           fieldNamesPar(iField)+" is missing.")
end
% Assert that all fields are function handles
if ~isa(ParamFunctions.fun_pLODinc,'function_handle')
    error("ParamFunctions.fun_pLODinc has to be a function handle.")
end
if ~isa(ParamFunctions.fun_pLODaz,'function_handle')
    error("ParamFunctions.fun_pLODaz has to be a function handle.")
end
if ~isa(ParamFunctions.fun_pLSD,'function_handle')
    error("ParamFunctions.fun_pLSD has to be a function handle.")
end

%% Target leaf area
% Assert that target leaf area is positive
if targetLeafArea <= 0
    error("targetLeafArea has to be positive.")
end
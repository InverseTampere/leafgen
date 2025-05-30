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

function successFlag = check_inputs_canopy_hull(TargetDistributions, ...
                                                LeafProperties, ...
                                                totalLeafArea)

%% Target distributions
% Assert the existence of all fields of the struct
fieldCheckTar = [isfield(TargetDistributions,'dTypeLADDh'), ...
                 isfield(TargetDistributions,'dTypeLADDd'), ...
                 isfield(TargetDistributions,'dTypeLADDc'), ...
                 isfield(TargetDistributions,'dTypeLODinc'), ...
                 isfield(TargetDistributions,'dTypeLODaz'), ...
                 isfield(TargetDistributions,'dTypeLSD'), ...
                 isfield(TargetDistributions,'pLADDh'), ...
                 isfield(TargetDistributions,'pLADDd'), ...
                 isfield(TargetDistributions,'pLADDc'), ...
                 isfield(TargetDistributions,'fun_pLODinc'), ...
                 isfield(TargetDistributions,'fun_pLODaz'), ...
                 isfield(TargetDistributions,'fun_pLSD'), ...
                 ];
fieldNamesTar = ["dTypeLADDh", ...
                 "dTypeLADDd", ...
                 "dTypeLADDc", ...
                 "dTypeLODinc", ...
                 "dTypeLODaz", ...
                 "dTypeLSD", ...
                 "pLADDh", ...
                 "pLADDd", ...
                 "pLADDc", ...
                 "fun_pLODinc", ...
                 "fun_pLODaz", ...
                 "fun_pLSD", ...
                  ];
for iField = 1:length(fieldCheckTar)
    assert(fieldCheckTar(iField),"TargetDistributions."...
           +fieldNamesTar(iField)+" is missing.")
end

%% Check the validity of LADD distribution names and parameters
% Relative height
dTypeH  = TargetDistributions.dTypeLADDh;
pLADDh = TargetDistributions.pLADDh;
if ~any(strcmp(dTypeH,{'uniform','polynomial','polynomialmixture', ...
        'weibull','weibullmixture','beta','betamixture',''}))
    error("LADD height distribution type not recognized.")
end
switch dTypeH
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'polynomial'
        % Assure that the polynomial gets only nonnegative values
        if any(polyval(pLADDh,0:0.001:1) < 0)
            error("TargetDistributions relative height polynomial gets"...
                  +" negative values on the interval [0,1].")
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
            error("TargetDistributions.pLADDh mixture model weight is"...
                  +" not on the interval [0,1].")
        end
        % Assure that the polynomial gets only nonnegative values
        if any(mmPolyValues < 0)
            error("TargetDistributions relative height polynomial gets"...
                  +" negative values on the interval [0,1].")
        end
    case 'weibull'
        % Assure that both parameters are positive
        if any(pLADDh <= 0)
            error("TargetDistributions.pLADDh all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
    case 'weibullmixture'
        % Assure that all Weibull parameters are positive
        if any(pLADDh(1:4) <= 0)
            error("TargetDistributions.pLADDh all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
        w = pLADDh(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.pLADDh mixture model weight is"...
                  +" not on the interval [0,1].")
        end
    case 'beta'
        % Assure that both parameters are positive
        if any(pLADDh <= 0)
            error("TargetDistributions.pLADDh all elements have to be"...
                  +" positive for beta distribution.")
        end
    case 'betamixture'
        % Assure that all beta parameters are positive
        if any(pLADDh(1:4) <= 0)
            error("TargetDistributions.pLADDh all elements have to be"...
                  +" positive for beta distribution.")
        end
        w = pLADDh(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.pLADDh mixture model weight is"...
                  +" not on the interval [0,1].")
        end
end
% Relative distance along subbranch
dTypeD  = TargetDistributions.dTypeLADDd;
pLADDd = TargetDistributions.pLADDd;
if ~any(strcmp(dTypeD,{'uniform','polynomial','polynomialmixture', ...
        'weibull','weibullmixture','beta','betamixture',''}))
    error("LADD distance from stem distribution type not recognized.")
end
switch dTypeD
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'polynomial'
        % Assure that the polynomial gets only nonnegative values
        if any(polyval(pLADDd,0:0.001:1) < 0)
            error("TargetDistributions relative distance along"...
                  +" subbranch polynomial gets negative values on the"...
                  +" interval [0,1].")
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
            error("TargetDistributions.pLADDd mixture model weight is"...
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
        if any(pLADDd <= 0)
            error("TargetDistributions.pLADDd all elements have to be"...
                  +" positive for truncated Weibull distribution.")
        end
    case 'weibullmixture'
        % Assure that all Weibull parameters are positive
        if any(pLADDd(1:4) <= 0)
            error("TargetDistributions.pLADDd all elements have to be"...
                  +" positive for truncated Weibull distributions.")
        end
        w = pLADDd(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.pLADDd mixture model weight is"...
                  +" not on the interval [0,1].")
        end
    case 'beta'
        % Assure that both parameters are positive
        if any(pLADDd <= 0)
            error("TargetDistributions.pLADDd all elements have to be"...
                  +" positive for beta distribution.")
        end
    case 'betamixture'
        % Assure that all beta parameters are positive
        if any(pLADDd(1:4) <= 0)
            error("TargetDistributions.pLADDd all elements have to be"...
                  +" positive for beta distributions.")
        end
        w = pLADDd(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.pLADDd mixture model weight is"...
                  +" not on the interval [0,1].")
        end
end
% Compass direction
dTypeC  = TargetDistributions.dTypeLADDc;
pLADDc = TargetDistributions.pLADDc;
if ~any(strcmp(dTypeC,{'uniform','vonmises','vonmisesmixture',''}))
    error("LADD compass direction distribution type not recognized.")
end
switch dTypeC
    case 'uniform'
        % parameters have no effect for uniform distribution
    case 'vonmises'
        % Assure that the second parameter is positive
        if pLADDc(2) <= 0
            error("TargetDistributions.pLADDc second parameter has"...
                  +" to be positive for von Mises distribution.")
        end
    case 'vonmisesmixturemodel'
        % Assure that second parameters are positive for both disrtibutions
        if pLADDc(2) <= 0 || pLADDc(4) <= 0
            error("TargetDistributions.pLADDc second and fourth"...
                  +" element have to be positive for von Mises"...
                  +" distributions.")
        end
        w = pLADDc(5); % mixture model weight
        % Check that the mixture model weight is between 0 and 1
        if w < 0 || w > 1
            error("TargetDistributions.pLADDc mixture model weight"...
                  +" is not on the interval [0,1].")
        end
end

%% Check validity of LOD parameter functions
% Assert that all fields are function handles
if ~isa(TargetDistributions.fun_pLODinc,'function_handle')
    error("TargetDistributions.fun_pLODinc has to be a function"...
          +" handle.")
end
if ~isa(TargetDistributions.fun_pLODaz,'function_handle')
    error("TargetDistributions.fun_pLODaz has to be a function handle.")
end

%% Check validity of LSD parameter function
if ~isa(TargetDistributions.fun_pLSD,'function_handle')
    error("TargetDistributions.fun_pLSD has to be a function"...
          +" handle.")
end

%% Check the leaf base geometry
if size(LeafProperties.vertices,2) ~= 3
    error("LeafProperties.vertices should have a size (nVertices x 3).")
end
if size(LeafProperties.triangles,2) ~= 3
    error("LeafProperties.triangles should have a size (nTriangles x 3).")
end
if max(LeafProperties.triangles,[],'all') > size(LeafProperties.vertices,1)
    error("LeafProperties.triangles refers to too large vertice indices.")
end

%% Total leaf area
% Assert that total leaf area is positive
if totalLeafArea <= 0
    error("totalLeafArea has to be positive.")
end

%% Mark success
successFlag = 1;

end
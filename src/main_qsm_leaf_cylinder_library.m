%% Clear variables and add paths

clear, close all
addpath("qsm-leaf-cylinder-library-method")
addpath("classes")
addpath("common-functions")
addpath("visualization")

%% Library LOD and LSD types

% Leaf orientation distribution types
LibraryDistributions.dTypeLODinc = 'dewit';
LibraryDistributions.dTypeLODaz  = 'vonmises';

% Leaf size distribution type
LibraryDistributions.dTypeLSD    = 'uniform';

%% Leaf distribution parameter nodes

% Inclination angle distribution nodes
Nodes.pLODinc1 = [-1 0 1];
Nodes.pLODinc2 = [2 3 4];

% Azimuth angle distribution nodes
Nodes.pLODaz1   = [0 pi/2 pi 3*pi/2];
Nodes.pLODaz2   = [0.01 0.5];

% Leaf size distribution nodes
Nodes.pLSD1 = [0.002 0.0025];
Nodes.pLSD2 = [0.003 0.0035];

%% Cylinder attribute nodes

Nodes.cylinderLength = [0.10 0.20];
Nodes.cylinderRadius = [0.01 0.05 0.1];
Nodes.cylinderInclinationAngle = [pi/4 3*pi/4];
Nodes.cylinderAzimuthAngle = [pi/4 3*pi/4 5*pi/4 7*pi/4];
Nodes.cylinderLeafArea = [0.2 0.5];

%% Leaf and petiole properties

% Vertices of the leaf basis geometry
LeafProperties.vertices = [-0.04  0.0   0.0;
                           0.0    0.08  0.0;
                           0.04   0.0   0.0];

% Triangles of the leaf basis geometry
LeafProperties.triangles = [1,  2,  3];

% Petiole length limits
LeafProperties.petioleLengthLimits = [0.05 0.1];

%% Generate leaf cylinder library
 
LeafCylinderLibrary = generate_leaf_cylinder_library( ...
                          LibraryDistributions,Nodes,LeafProperties);

%% Save the library for later use if needed

saveLibrary = false;

if saveLibrary == true
    save('NewLeafCyliderLibrary.mat','-struct','LeafCylinderLibrary')
end

%% Initialize QSM

filename = "example-data/ExampleQSM.mat";
QSM = importdata(filename);

%% Define target LADD marginal distributions

% LADD relative height
TargetLADD.dTypeLADDh = 'beta';
TargetLADD.hParams = [4.5 1];

% LADD relative distance along sub-branch
TargetLADD.dTypeLADDd = 'beta';
TargetLADD.dParams = [7 1];

% LADD compass direction
TargetLADD.dTypeLADDc = 'vonmises';
TargetLADD.cParams = [pi 0.1];

%% Define parameter functions for LOD and LSD 

% LOD inclination angle
ParamFunctions.fun_inc_params = @(h,d,c) [-1,4];

% LOD azimuth angle
ParamFunctions.fun_az_params = @(h,d,c) [3.3, 0.25];

% LSD
ParamFunctions.fun_size_params = @(h,d,c) [0.0021, 0.0038];

%% Populate QSM with leaves using leaf cylinder library

targetLeafArea = 50;

[Leaves,QSMbc] = transform_leaf_cylinders(QSM,LeafCylinderLibrary, ...
                                          TargetLADD,ParamFunctions, ...
                                          targetLeafArea);

%% Visualize the QSM with generated foliage

figure(1), clf
% Plot QSM
hQSM = QSMbc.plot_model();
% Set bark color
set(hQSM,'FaceColor',[150,100,50]./255,'EdgeColor','none');

hold on;

% Plot leaves
hLeaf = Leaves.plot_leaves();
% Set leaf color
set(hLeaf,'FaceColor',[0,150,0]./255,'EdgeColor','none');

hold off;
axis equal;
xlabel('x')
ylabel('y')
zlabel('z')

%% Plot LADD marginal distributions

plot_LADD_h_QSM(QSMbc,Leaves,TargetLADD);
plot_LADD_d_QSM(QSMbc,Leaves,TargetLADD);
plot_LADD_c_QSM(QSMbc,Leaves,TargetLADD);

%% Export leaves and QSM in OBJ-format

% Precision parameter for export
precision = 5;

% Exporting to obj files
Leaves.export_geometry('OBJ',true,'leaves_export.obj',precision);
QSMbc.export('OBJ','qsm_export.obj','Precision',precision);
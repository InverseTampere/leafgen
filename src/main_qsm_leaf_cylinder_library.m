%% Clear variables and add paths

clear, close all
addpath("qsm-leaf-cylinder-library-method")
addpath("classes")
addpath("common-functions")
addpath("visualization")

%% Leaf and petiole properties

% Vertices of the leaf base geometry
LeafProperties.vertices = [0.0    0.0    0.0;
                           -0.04  0.02   0.0;
                           0.0    0.10   0.0;
                           0.04   0.02   0.0];

% Triangles of the leaf base geometry
LeafProperties.triangles = [1 2 3;
                            1 3 4];

% Petiole length limits
LeafProperties.petioleLengthLimits = [0.08 0.10];

%% Library LOD and LSD types

% Leaf orientation distribution types
LibraryDistributions.dTypeLODinc = 'dewit';
LibraryDistributions.dTypeLODaz  = 'vonmises';

% Leaf size distribution type
LibraryDistributions.dTypeLSD    = 'uniform';

%% Leaf distribution parameter nodes

% Inclination angle distribution nodes
Nodes.pLODinc1 = -1%[-1 0 1];
Nodes.pLODinc2 = 4%[2 3 4];

% Azimuth angle distribution nodes
Nodes.pLODaz1 = pi%[0 pi/2 pi 3*pi/2];
Nodes.pLODaz2 = 0.01%[0.01 0.5];

% Leaf size distribution nodes
Nodes.pLSD1 = 0.0025%[0.002 0.0025];
Nodes.pLSD2 = 0.0030%[0.003 0.0035];

%% Cylinder attribute nodes

Nodes.cylinderLength = [0.10 0.20];
Nodes.cylinderRadius = [0.01 0.05];
Nodes.cylinderInclinationAngle = [pi/4 3*pi/4];
Nodes.cylinderAzimuthAngle = [pi/4 3*pi/4 5*pi/4 7*pi/4];
Nodes.cylinderLeafArea = [0.2 0.5 0.8];

%% Generate leaf cylinder library
 
LeafCylinderLibrary = generate_leaf_cylinder_library( ...
                          LibraryDistributions,Nodes,LeafProperties);

%% Save the library for later use

save("NewLeafCyliderLibrary.mat",'-struct','LeafCylinderLibrary')

%% Initialize QSM

filename = "example-data/ExampleQSM.mat";
QSM = importdata(filename);

%% Define target LADD marginal distributions

% LADD relative height
TargetLADD.dTypeLADDh = 'beta';
TargetLADD.pLADDh = [22 3];

% LADD relative distance along sub-branch
TargetLADD.dTypeLADDd = 'weibull';
TargetLADD.pLADDd = [3.3 2.8];

% LADD compass direction
TargetLADD.dTypeLADDc = 'vonmises';
TargetLADD.pLADDc = [5/4*pi 0.1];

%% Define parameter functions for LOD and LSD 

% LOD inclination angle
ParamFunctions.fun_pLODinc = @(h,d,c) [-1 4];

% LOD azimuth angle
ParamFunctions.fun_pLODaz = @(h,d,c) [pi 0.01];

% LSD
ParamFunctions.fun_pLSD = @(h,d,c) [0.0025, 0.0030];

%% Set the target leaf area

totalLeafArea = 50;

%% Populate QSM with leaves using leaf cylinder library

[Leaves,QSMbc] = transform_leaf_cylinders(QSM,LeafCylinderLibrary, ...
                                          TargetLADD,ParamFunctions, ...
                                          totalLeafArea);

%% Visualize the QSM with generated foliage

% Initialize figure
figure(1), clf, hold on

% Plot leaves
hLeaf = Leaves.plot_leaves();

% Set leaf color
set(hLeaf,'FaceColor',[0 150 0]./255,'EdgeColor','none');

% Plot QSM
hQSM = QSMbc.plot_model();

% Set bark color
set(hQSM,'FaceColor',[150 100 50]./255,'EdgeColor','none');

% Set figure properties
hold off;
axis equal;
xlabel('x')
ylabel('y')
zlabel('z')

%% Plot LADD marginal distributions

plot_LADD_h_QSM(QSMbc,Leaves,TargetLADD);
plot_LADD_d_QSM(QSMbc,Leaves,TargetLADD);
plot_LADD_c_QSM(QSMbc,Leaves,TargetLADD);

%% Plot LOD marginal distributions

plot_LOD_inc_QSM(QSMbc,Leaves);
plot_LOD_az_QSM(QSMbc,Leaves);

%% Plot LSD

plot_LSD_QSM(QSMbc,Leaves);

%% Export leaves and QSM in OBJ-format

% Precision parameter for export
precision = 5;

% Exporting to obj files
Leaves.export_geometry('OBJ',true,'leaves_export.obj',precision);
QSMbc.export('OBJ','qsm_export.obj','Precision',precision);
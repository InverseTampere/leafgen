%% Clear variables and add necessary paths
clear, close all
addpath('classes/');

%% Leaf distribution types

% Leaf orientation distribution types
LibraryDistributions.dTypeLODinc = 'dewit';
LibraryDistributions.dTypeLODaz  = 'uniform'; %'vonmises';

% Leaf size distribution type
LibraryDistributions.dTypeLSD    = 'uniform';

%% Leaf distribution parameter nodes

% Inclination angle distribution nodes
nIncNodes1    = 1;
nIncNodes2    = 1;
aInterval     = [-1 -1]; %[-1 1];
bInterval     = [ 4  4]; %[ 2 4];
Nodes.pLODinc1 = linspace(aInterval(1),aInterval(2),nIncNodes1);
Nodes.pLODinc2 = linspace(bInterval(1),bInterval(2),nIncNodes2);

% Azimuth angle distribution nodes
% nAzNodes1      = 3;
% nAzNodes2      = 2;
% muInterval     = [0 2*pi]; % *(1-1/nAzNodes1)
% kappaInterval  = [0.01 0.5];
% Nodes.pLODaz1   = linspace(muInterval(1),muInterval(2),nAzNodes1);
% Nodes.pLODaz2   = linspace(kappaInterval(1),kappaInterval(2),nAzNodes2);

% Leaf size distribution nodes
nLSDNodes1 = 1;
nLSDNodes2 = 1;
lbInterval = [0.002 0.002];
ubInterval = [0.004 0.004];
Nodes.pLSD1 = linspace(lbInterval(1),lbInterval(2),nLSDNodes1);
Nodes.pLSD2 = linspace(ubInterval(1),ubInterval(2),nLSDNodes2);

%% Cylinder attribute nodes

nCylLenNodes = 3;
nCylRadNodes = 3;
nCylIncNodes = 2;
nCylAzNodes  = 4;
nCylArNodes  = 2;

Nodes.cylinderLength = linspace(0.10,0.20,nCylLenNodes);
Nodes.cylinderRadius = linspace(0.01,0.10,nCylRadNodes);
Nodes.cylinderInclinationAngle = linspace(0.5*pi/nCylIncNodes, ...
                                          pi-0.5*pi/nCylIncNodes, ...
                                          nCylIncNodes);
Nodes.cylinderAzimuthAngle = linspace(0.5*2*pi*nCylAzNodes, ...
                                      2*pi-0.5*2*pi*nCylAzNodes, ...
                                      nCylAzNodes);
Nodes.cylinderLeafArea = linspace(0.2,0.5,nCylArNodes);

%% Leaf and twig properties

% Twig length limits
LeafProperties.twigLengthLimits = [0.05 0.1];

% Vertices of the leaf basis geometry
LeafProperties.vertices = [-0.04  0.0   0.0;
                           0.0    0.08  0.0;
                           0.04   0.0   0.0];

% Triangles of the leaf basis geometry
LeafProperties.triangles = [1,  2,  3];

%% Generate leaf-cylinder library
% LeafCylinderLibrary = load('LeafCyliderLibraryExample.mat');

tic
LeafCylinderLibrary = generate_leaf_cylinder_library(Nodes, ...
                        LibraryDistributions, ...
                        LeafProperties, ...
                        'nLeafObjectsPerNode',10, ...
                        'PreventIntersections',false);
toc

%% Save the leaf cylinder library as a .mat file

% save('NewLeafCyliderLibraryExample.mat','-struct','LeafCylinderLibrary')
% disp('---------------------')
% disp(LeafCylinderLibrary.totalNodes)
% s = dir('NewLeafCyliderLibraryExample.mat');         
% disp(s.bytes)
% return

%% Initialize QSM object.
% QSM = QSMBCylindrical('example');
filename = "qsm_Small.mat";
load(filename); % contains the struct named "qsm"
QSM = QSMBCylindrical(qsm);

%% Define target leaf distributions

% LADD relative height
TargetLADD.dTypeLADDh = 'beta';
TargetLADD.hParams = [4.5 1];

% LADD relative distance along sub-branch
TargetLADD.dTypeLADDd = 'beta';
TargetLADD.dParams = [7 1];

% LADD compass direction
TargetLADD.dTypeLADDc = 'vonmises';
TargetLADD.cParams = [pi 0.1];

% LOD inclination angle
ParamFunctions.fun_inc_params = @(h,d,c) [-1,4];

% LOD azimuth angle
ParamFunctions.fun_az_params = @(h,d,c) [3.3, 0.25];

% LSD
ParamFunctions.fun_size_params = @(h,d,c) [0.0021, 0.0038];

%% Populate QSM with leaves using leaf-clinder library

% Target area for the leaves
targetLeafArea = 50;

tic
Leaves = populate_qsm_with_leaves(QSM,LeafCylinderLibrary, ...
                                  TargetLADD,ParamFunctions, ...
                                  targetLeafArea);
toc

%% Visualize the QSM with generated foliage

figure(1), clf
% Plot QSM
hQSM = QSM.plot_model();
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

%% Plot leaf distributions

plot_LADD_h_LCL(QSM,Leaves,TargetLADD);
plot_LADD_d_LCL(QSM,Leaves,TargetLADD);
plot_LADD_c_LCL(QSM,Leaves,TargetLADD);

%% Plots for debugging LOD and LSD

plot_LOD_inc_LCL(QSM,Leaves, ...
                 LeafCylinderLibrary.LeafDistributions, ...
                 ParamFunctions,10);
% plot_LOD_az_LCL
% plot_LSD_LCL





%% Clear variables and add necessary paths
clear, close all
addpath('classes/');

%% Leaf distribution types

% Leaf orientation distribution types
LibraryDistributions.dTypeLOD_inc = 'dewit';
LibraryDistributions.dTypeLOD_az  = 'uniform'; %'vonmises';

% Leaf size distribution type
LibraryDistributions.dTypeLSD    = 'uniform';

%% Leaf distribution parameter nodes

% Inclination angle distribution nodes
nIncNodes1    = 1;
nIncNodes2    = 1;
aInterval     = [-1 -1]; %[-1 1];
bInterval     = [ 4  4]; %[ 2 4];
Nodes.LodInc1 = linspace(aInterval(1),aInterval(2),nIncNodes1);
Nodes.LodInc2 = linspace(bInterval(1),bInterval(2),nIncNodes2);

% Azimuth angle distribution nodes
% nAzNodes1      = 3;
% nAzNodes2      = 2;
% muInterval     = [0 2*pi]; % *(1-1/nAzNodes1)
% kappaInterval  = [0.01 0.5];
% Nodes.LodAz1   = linspace(muInterval(1),muInterval(2),nAzNodes1);
% Nodes.LodAz2   = linspace(kappaInterval(1),kappaInterval(2),nAzNodes2);

% Leaf size distribution nodes
nLsdNodes1 = 1;
nLsdNodes2 = 1;
lbInterval = [0.002 0.002];
ubInterval = [0.004 0.004];
Nodes.Lsd1 = linspace(lbInterval(1),lbInterval(2),nLsdNodes1);
Nodes.Lsd2 = linspace(ubInterval(1),ubInterval(2),nLsdNodes2);

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

%% Leaf and twig base parameters

twigLengthLimits = [0.05 0.1];

% Vertices of the leaf basis geometry
vertices = [
    -0.04  0.0   0.0;
    0.0    0.08  0.0;
    0.04   0.0   0.0
];

% Triangles of the leaf basis geometry
tris = [
     1,  2,  3
];

%% Generate leaf-cylinder library
% LeafCylinderLibrary = load('LeafCyliderLibraryExample.mat');

tic
LeafCylinderLibrary = generate_leaf_cylinder_library(Nodes, ...
                        LibraryDistributions, ...
                        twigLengthLimits, ...
                        vertices, ...
                        tris, ...
                        'nLeafObjectsPerNode',1, ...
                        'PreventIntersections',true);
toc

%% Initialize QSM object.
% QSM = QSMBCylindrical('example');
filename = "qsm_Small.mat";
load(filename); % contains the struct named "qsm"
QSM = QSMBCylindrical(qsm);

%% Define target leaf distributions

% LADD relative height
TargetLADD.dTypeLADD_h = 'beta';
TargetLADD.p_h = [4.5 1];
TargetLADD.nBins_h = 10;

% LADD relative distance along sub-branch
TargetLADD.dTypeLADD_d = 'beta';
TargetLADD.p_d = [7 1];
TargetLADD.nBins_d = 10;

% LADD compass direction
TargetLADD.dTypeLADD_c = 'vonmises';
TargetLADD.p_c = [pi 0.1];
TargetLADD.nBins_c = 10;

% LOD inclination angle
ParamFunctions.fun_inc_params = @(h,d,c) [1,2];

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
set(hLeaf,'FaceColor',[120,150,80]./255,'EdgeColor','none');

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
                 ParamFunctions,10)
% plot_LOD_az_LCL
% plot_LSD_LCL

%% Save the leaf cylinder library as a .mat file

return
save('NewLeafCyliderLibraryExample.mat','-struct','LeafCylinderLibrary')




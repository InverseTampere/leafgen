%% Clear variables and add necessary paths
clear, close all
addpath('classes/');

%% Leaf distribution types

% Leaf orientation distribution types
LeafDistributions.dTypeLodInc = 'dewit';
LeafDistributions.dTypeLodAz  = 'vonmises';

% Leaf size distribution types
LeafDistributions.dTypeLsd    = 'uniform';

%% Leaf-cylinder library nodes

% Inclination angle distribution nodes
nIncNodes1    = 2;
nIncNodes2    = 2;
aInterval     = [-1 1];
bInterval     = [ 2 4];
Nodes.LodInc1 = linspace(aInterval(1),aInterval(2),nIncNodes1);
Nodes.LodInc2 = linspace(bInterval(1),bInterval(2),nIncNodes2);

% Azimuth angle distribution nodes
nAzNodes1      = 3;
nAzNodes2      = 2;
muInterval     = [0 2*pi]; % *(1-1/nAzNodes1)
kappaInterval  = [0.01 0.5];
Nodes.LodAz1   = linspace(muInterval(1),muInterval(2),nAzNodes1);
Nodes.LodAz2   = linspace(kappaInterval(1),kappaInterval(2),nAzNodes2);

% Leaf size distribution nodes
nLsdNodes1 = 2;
nLsdNodes2 = 2;
lbInterval = [0.002 0.0025];
ubInterval = [0.0035 0.004];
Nodes.Lsd1 = linspace(lbInterval(1),lbInterval(2),nLsdNodes1);
Nodes.Lsd2 = linspace(ubInterval(1),ubInterval(2),nLsdNodes2);

%% Cylinder attribute nodes

nCylLenNodes = 2;
nCylRadNodes = 2;
nCylIncNodes = 3;
nCylAzNodes  = 3;
nCylArNodes  = 2;

Nodes.cylinderLength = linspace(0.10,0.5,nCylLenNodes);
Nodes.cylinderRadius = linspace(0.05,0.2,nCylRadNodes);
Nodes.cylinderInclinationAngle = linspace(0,pi,nCylIncNodes);
Nodes.cylinderAzimuthAngle = linspace(0,2*pi*(1-(1/nCylAzNodes)),nCylAzNodes);
Nodes.cylinderLeafArea = linspace(0.3,1,nCylArNodes);

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
LeafCylinderLibrary = load('LeafCyliderLibraryExample.mat');

% tic
% LeafCylinderLibrary = generate_leaf_cylinder_library(Nodes, ...
%                         LeafDistributions, ...
%                         twigLengthLimits, ...
%                         vertices, ...
%                         tris, ...
%                         'nLeafObjectsPerNode',1, ...
%                         'PreventIntersections',true);
% toc

%% Initialize QSM object.
% QSM = QSMBCylindrical('example');
filename = "qsm_Small.mat";
load(filename); % contains the struct named "qsm"
QSM = QSMBCylindrical(qsm);

%% Define target leaf distributions

% LADD relative height
TargetDistributions.dType_h = 'beta';
TargetDistributions.p_h = [4.5 1];
TargetDistributions.nBins_h = 10;

% LADD relative distance along sub-branch
TargetDistributions.dType_d = 'beta';
TargetDistributions.p_d = [7 1];
TargetDistributions.nBins_d = 10;

% LADD compass direction
TargetDistributions.dType_c = 'vonmises';
TargetDistributions.p_c = [pi 0.1];
TargetDistributions.nBins_c = 10;

% LOD inclination angle
TargetDistributions.dType_inc = 'dewit';
TargetDistributions.fun_inc_params = @(h,d,c) [1,2];
% TargetDistributions.p_inc = [1 2];

% LOD azimuth angle
TargetDistributions.dType_az = 'vonmises';
TargetDistributions.fun_az_params = @(h,d,c) [3.3, 0.25];
% TargetDistributions.p_az = [3.3 0.25];

% LSD
TargetDistributions.dType_size = 'uniform';
TargetDistributions.fun_size_params = @(h,d,c) [0.0021, 0.0038];
% TargetDistributions.p_size =  [0.0021 0.0038];

%% Populate QSM with leaves using leaf-clinder library

% Target area for the leaves
targetLeafArea = 50;

tic
Leaves = populate_qsm_with_leaves(QSM,LeafCylinderLibrary, ...
                                  TargetDistributions,targetLeafArea);
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

plot_LADD_h_LCL(QSM,Leaves,TargetDistributions);
plot_LADD_d_LCL(QSM,Leaves,TargetDistributions);
plot_LADD_c_LCL(QSM,Leaves,TargetDistributions);

%% Save the leaf cylinder library as a .mat file

return
save('BigLeafCyliderLibraryExample.mat','-struct','LeafCylinderLibrary')




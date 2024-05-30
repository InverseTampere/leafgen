clear, close all
addpath("C:\Users\bppimo\Documents\GitHub\foliage-generation-repo\src")
addpath("C:\Users\bppimo\Documents\GitHub\foliage-generation-repo\src\classes")

%% Initialize leaf base geometry

% Vertices of the leaf basis geometry
LeafProperties.vertices = [-0.04  0.0   0.0;
                           0.0    0.08  0.0;
                           0.04   0.0   0.0];

% Triangles of the leaf basis geometry
LeafProperties.triangles = [1, 2, 3];

LeafProperties.twigLengthLimits = [0.05 0.1];

%% Define target leaf distributions

% LADD relative height
TargetDistributions.dTypeLADDh = 'betamixturemodel';
TargetDistributions.hParams = [22 3 41 50 0.85];

% LADD relative distance from stem
TargetDistributions.dTypeLADDd = 'beta';
TargetDistributions.dParams = [2 1];

% LADD compass direction
TargetDistributions.dTypeLADDc = 'vonmisesmixturemodel';
TargetDistributions.cParams = [pi 0.1 6/5*pi 0.1 0.6];

% LOD inclination angle
TargetDistributions.dTypeLODinc = 'dewit';
TargetDistributions.fun_inc_params = @(h,d,c) [1,2];

% LOD azimuth angle
TargetDistributions.dTypeLODaz = 'vonmises';
TargetDistributions.fun_az_params = @(h,d,c) [3.3, 0.25];

% LSD
TargetDistributions.dTypeLSD = 'uniform';
TargetDistributions.fun_size_params = @(h,d,c) [0.0021, 0.0038];

%% Initialize QSM object.
% QSM = QSMBCylindrical('example');
filename = "qsm_Small.mat";
load(filename); % contains the struct named "qsm"
QSM = QSMBCylindrical(qsm);

%% Generate foliage on QSM
totalLeafArea = 50;
tic
Leaves = generate_foliage_qsm_direct(QSM,TargetDistributions, ...
                                     LeafProperties,totalLeafArea);
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

plot_LADD_h_LCL(QSM,Leaves,TargetDistributions);
plot_LADD_d_LCL(QSM,Leaves,TargetDistributions);
plot_LADD_c_LCL(QSM,Leaves,TargetDistributions);

%% Plots for debugging LOD and LSD

plot_LOD_inc_LCL(QSM,Leaves,TargetDistributions,TargetDistributions,10);
% plot_LOD_az_LCL
% plot_LSD_LCL
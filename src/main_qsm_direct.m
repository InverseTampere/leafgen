%% Clear variables and add paths

clear, close all
addpath("qsm direct method")
addpath("classes")
addpath("common functions")
addpath("visualization")

%% Initialize QSM

filename = "example data/ExampleQSM.mat";
qsm = importdata(filename);
QSM = QSMBCylindrical(qsm);

%% Initialize leaf base geometry

% Vertices of the leaf basis geometry
LeafProperties.vertices = [-0.04  0.0   0.0;
                           0.0    0.08  0.0;
                           0.04   0.0   0.0];

% Triangles of the leaf basis geometry
LeafProperties.triangles = [1, 2, 3];

% Twig length limits
LeafProperties.twigLengthLimits = [0.05 0.1];

%% Define target leaf distributions

% LADD relative height
TargetDistributions.dTypeLADDh = 'betamixture';
TargetDistributions.hParams = [22 3 41 50 0.85];

% LADD relative distance from stem
TargetDistributions.dTypeLADDd = 'beta';
TargetDistributions.dParams = [2 1];

% LADD compass direction
TargetDistributions.dTypeLADDc = 'vonmises';
TargetDistributions.cParams = [5/4*pi 0.1];

% LOD inclination angle
TargetDistributions.dTypeLODinc = 'dewit';
TargetDistributions.fun_inc_params = @(h,d,c) [1,2];

% LOD azimuth angle
TargetDistributions.dTypeLODaz = 'uniform';
TargetDistributions.fun_az_params = @(h,d,c) [];

% LSD
TargetDistributions.dTypeLSD = 'uniform';
TargetDistributions.fun_size_params = @(h,d,c) [0.0021, 0.0038];

%% Generate foliage on QSM

totalLeafArea = 50;

Leaves = generate_foliage_qsm_direct(QSM,TargetDistributions, ...
                                     LeafProperties,totalLeafArea);

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

%% Plot LADD marginal distributions

plot_LADD_h_QSM(QSM,Leaves,TargetDistributions);
plot_LADD_d_QSM(QSM,Leaves,TargetDistributions);
plot_LADD_c_QSM(QSM,Leaves,TargetDistributions);

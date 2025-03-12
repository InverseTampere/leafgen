%% Clear variables and add paths

clear, close all
addpath("qsm-direct-method")
addpath("classes")
addpath("common-functions")
addpath("visualization")

%% Initialize QSM

filename = "example-data/ExampleQSM.mat";
QSM = importdata(filename);

%% Initialize leaf base geometry

% Vertices of the leaf base geometry
LeafProperties.vertices = [0.0    0.0    0.0;
                           -0.04  0.02   0.0;
                           0.0    0.10   0.0;
                           0.04   0.02   0.0];

% Triangles of the leaf base geometry
LeafProperties.triangles = [1 2 3;
                            1 3 4];

%% Define petiole length sampling interval

LeafProperties.petioleLengthLimits = [0.08 0.10];

%% Define target leaf distributions

% LADD relative height
TargetDistributions.dTypeLADDh = 'beta';
TargetDistributions.pLADDh = [22 3];

% LADD relative branch distance
TargetDistributions.dTypeLADDd = 'weibull';
TargetDistributions.pLADDd = [3.3 2.8];

% LADD compass direction
TargetDistributions.dTypeLADDc = 'vonmises';
TargetDistributions.pLADDc = [5/4*pi 0.1];

% LOD inclination angle
TargetDistributions.dTypeLODinc = 'dewit';
TargetDistributions.fun_inc_params = @(h,d,c) [1 2];

% LOD azimuth angle
TargetDistributions.dTypeLODaz = 'uniform';
TargetDistributions.fun_az_params = @(h,d,c) [];

% LSD
TargetDistributions.dTypeLSD = 'uniform';
TargetDistributions.fun_size_params = @(h,d,c) [0.0021 0.0038];

%% Generate foliage on QSM

totalLeafArea = 50;

[Leaves,QSMbc] = generate_foliage_qsm_direct(QSM,TargetDistributions, ...
                                             LeafProperties,totalLeafArea);

%% Visualize the QSM with generated foliage

figure(1), clf
% Plot QSM
hQSM = QSMbc.plot_model();
% Set bark color
set(hQSM,'FaceColor',[150 100 50]./255,'EdgeColor','none');

hold on;

% Plot leaves
hLeaf = Leaves.plot_leaves();
% Set leaf color
set(hLeaf,'FaceColor',[0 150 0]./255,'EdgeColor','none');

hold off;
axis equal;
xlabel('x')
ylabel('y')
zlabel('z')

%% Plot LADD marginal distributions

plot_LADD_h_QSM(QSMbc,Leaves,TargetDistributions);
plot_LADD_d_QSM(QSMbc,Leaves,TargetDistributions);
plot_LADD_c_QSM(QSMbc,Leaves,TargetDistributions);

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
%% Clear variables and add paths

clear, close all
addpath("canopy-hull-method")
addpath("classes")
addpath("common-functions")
addpath("visualization")

%% Initialize point cloud

filename = 'example-data/examplePC.mat';
ptCloud = importdata(filename);
ptCloud = double(ptCloud);

%% Initialize leaf base geometry

% Vertices of the leaf base geometry
LeafProperties.vertices = [0.0    0.0    0.0;
                           -0.04  0.02   0.0;
                           0.0    0.10   0.0;
                           0.04   0.02   0.0];

% Triangles of the leaf base geometry
LeafProperties.triangles = [1 2 3;
                            1 3 4];

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
TargetDistributions.fun_pLODinc = @(h,d,c) [1 2];

% LOD azimuth angle
TargetDistributions.dTypeLODaz = 'uniform';
TargetDistributions.fun_pLODaz = @(h,d,c) [];

% LSD
TargetDistributions.dTypeLSD = 'normal';
TargetDistributions.fun_pLSD = @(h,d,c) [0.004 0.00025^2];


%% Define stem location

% Translate the lowest point of the cloud to the plane z=0
ptCloud = ptCloud - min(ptCloud(:,3));

% Translate the trunk of the tree to origin
tfBottom = ptCloud(:,3) < 1;
ptCloud = ptCloud - [mean(ptCloud(tfBottom,1:2)) 0];

% Set the coordinates for the stem
stemCoordinates = [   0        0                   0;
                   0.47   -0.125                11.5;
                      0        0   max(ptCloud(:,3))];

%% Set the target leaf area

totalLeafArea = 25;

%% Generate foliage inside point cloud

[Leaves,aShape] = generate_foliage_canopy_hull(ptCloud, ...
                                     TargetDistributions, ...
                                     LeafProperties,totalLeafArea, ...
                                     'StemCoordinates',stemCoordinates);

%% Visualize the foliage

% Initialize figure
figure, clf
tiledlayout(1,2)

% Initialize first subplot
ax1 = nexttile;

% Plot point cloud
pc = aShape.Points;
plot3(pc(:,1),pc(:,2),pc(:,3),'k.','MarkerSize',1)

% Set subfigure properties
hold on, grid on, axis equal
xlabel('x')
ylabel('y')
zlabel('z')

% Plot alpha shape
plot(aShape,'FaceColor','m','FaceAlpha',0.2)

% Plot the defined stem
if exist('stemCoordinates','var')
    plot3(stemCoordinates(:,1),stemCoordinates(:,2), ...
        stemCoordinates(:,3),'c-','LineWidth',3);
else
    pcTop = max(pc(:,3));
    plot3([0 0],[0 0],[0 pcTop],'c-','LineWidth',3)
end

% Save the coordinate limits
xl = xlim; 
yl = ylim; 
zl = zlim;

% Initialize second subplot
ax2 = nexttile;

% Plot leaves
hLeaf = Leaves.plot_leaves();

% Set leaf color
set(hLeaf,'FaceColor',[0,150,0]./255,'EdgeColor','none');

% Set subfigure properties
grid on, axis equal, xlim(ax2,xl), ylim(ax2,yl), zlim(ax2,zl)
xlabel('x')
ylabel('y')
zlabel('z')

% Link the point of view for the subfigures
Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', ...
    'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

%% Plot LADD marginal distributions

plot_LADD_h_CH(aShape,Leaves,TargetDistributions);
plot_LADD_d_CH(aShape,Leaves,TargetDistributions, ...
               'StemCoordinates',stemCoordinates);
plot_LADD_c_CH(aShape,Leaves,TargetDistributions, ...
               'StemCoordinates',stemCoordinates);

%% Plot LOD marginal distributions

plot_LOD_inc_CH(aShape,Leaves,'StemCoordinates',stemCoordinates);
plot_LOD_az_CH(aShape,Leaves,'StemCoordinates',stemCoordinates);

%% Plot LSD

plot_LSD_CH(aShape,Leaves,'StemCoordinates',stemCoordinates);

%% Export leaves in OBJ-format

% Precision parameter for export
precision = 5;

% Exporting to obj file
Leaves.export_geometry('OBJ',true,'leaves_export.obj',precision);
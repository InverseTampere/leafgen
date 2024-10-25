%% Clear variables and add paths

clear, close all
addpath("canopy hull method")
addpath("classes")
addpath("common functions")
addpath("visualization")

%% Point cloud

filename = 'example data/examplePC.mat';
pCloud = importdata(filename);
pCloud = double(pCloud);

% Translate the lowest point of the cloud to the plane z=0
pCloud = pCloud - min(pCloud(:,3));

% Translate the trunk of the tree to origin
tfBottom = pCloud(:,3) < 1;
pCloud = pCloud - [mean(pCloud(tfBottom,1:2)) 0];

%% Initialize leaf base geometry

% Vertices of the leaf basis geometry
LeafProperties.vertices = [-0.04  0.0   0.0;
                           0.0    0.08  0.0;
                           0.04   0.0   0.0];

% Triangles of the leaf basis geometry
LeafProperties.triangles = [1, 2, 3];

%% Define target leaf distributions

% LADD relative height
TargetDistributions.dTypeLADDh = 'betamixture';
TargetDistributions.hParams = [22 3 41 50 0.85];

% LADD relative distance from stem
TargetDistributions.dTypeLADDd = 'beta';
TargetDistributions.dParams = [2 1];

% LADD compass direction
TargetDistributions.dTypeLADDc = 'vonmisesmixture';
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


%% Stem coordinates

stemCoordinates = [   0,      0,                  0;
                   0.47, -0.125,               11.5;
                      0,      0,   max(pCloud(:,3))];

%% Generate foliage

totalLeafArea = 25;

[Leaves,aShape] = generate_foliage_canopy_hull(pCloud, ...
                                     TargetDistributions, ...
                                     LeafProperties,totalLeafArea, ...
                                     'StemCoordinates',stemCoordinates);

%% Visualize the foliage

figure, clf

% Plot point cloud, alphashape and stem
tiledlayout(1,2)
ax1 = nexttile;
pc = aShape.Points;
plot3(pc(:,1),pc(:,2),pc(:,3),'k.','MarkerSize',1)
hold on, grid on, axis equal
plot(aShape,'FaceColor','m','FaceAlpha',0.2)
pcTop = max(pc(:,3));
if exist('stemCoordinates','var')
    plot3(stemCoordinates(:,1),stemCoordinates(:,2), ...
        stemCoordinates(:,3),'c-','LineWidth',3);
else
    plot3([0 0],[0 0],[0 pcTop],'c-','LineWidth',3)
end
xl = xlim;
yl = ylim;
zl = zlim;
ax2 = nexttile;
Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', ...
    'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

% Plot leaves
hLeaf = Leaves.plot_leaves();
% Set leaf color
set(hLeaf,'FaceColor',[0,150,0]./255,'EdgeColor','none');
grid on, axis equal, xlim(ax2,xl), ylim(ax2,yl), zlim(ax2,zl)
xlabel('x')
ylabel('y')
zlabel('z')

%% Plot LADD marginal distributions

plot_LADD_h_CH(aShape,Leaves,TargetDistributions);
plot_LADD_d_CH(aShape,Leaves,TargetDistributions, ...
               'StemCoordinates',stemCoordinates);
plot_LADD_c_CH(aShape,Leaves,TargetDistributions, ...
               'StemCoordinates',stemCoordinates);

%% Export leaves in OBJ-format

% Precision parameter for export
precision = 5;

% Exporting to obj file
Leaves.export_geometry('OBJ',false,'leaves_exportF.obj',precision);
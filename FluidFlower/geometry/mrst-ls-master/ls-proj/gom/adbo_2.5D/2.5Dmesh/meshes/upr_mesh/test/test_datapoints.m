%% Make datapoints for mesh conforming to boundary, fault, folded layers
%
% rock.perm(ucids.unit_cell_ids{33},1) = 700*(milli*darcy);
% plotToolbar(G, rock.perm(:,1)/(milli*darcy)); axis equal off; view([90 0])
% colormap(turbo)
% ylim([11000 12000]); zlim([0 250])
%

% Boundary
pth = fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/2.5Dmesh/meshes/upr_mesh');
h = 5000; % height (m)
l = 20000; % length
stratiPoints = struct('boundary', [], 'lines', [], 'wells', []);
stratiPoints.boundary =      [0, 0;
                              l, 0;
                              l, h; ...
                              0, h; ...
                              0, 0];
                          
% Wells (to have cells with correct centroids)
stratiPoints.wells{1} = [10000, 2375];              %I1 (HW)

% Lines (1 fault, 2 lines are shear zone boundaries)
stratiPoints.lines{1} = [0 1000;
                         20000 1000];
stratiPoints.lines{2} = [0 1500;
                         20000 1500];
                     
% subdivisions
stratiPoints.lines{3} = [0     1050;
                         20000 1050];
stratiPoints.lines{4} = [0     1100;
                         20000 1100];
stratiPoints.lines{5} = [0     1150;
                         20000 1150];
stratiPoints.lines{6} = [0     1200;
                         20000 1200];
                     
% Save data
save('test_datapoints.mat', 'stratiPoints');
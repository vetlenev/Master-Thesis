%% Make datapoints for mesh conforming to boundary, fault, folded layers
%
% rock.perm(ucids.unit_cell_ids{33},1) = 700*(milli*darcy);
% plotToolbar(G, rock.perm(:,1)/(milli*darcy)); axis equal off; view([90 0])
% colormap(turbo)
% ylim([11000 12000]); zlim([0 250])
%
clear, close all

% Boundary
pth = fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/2.5Dmesh/meshes/upr_mesh/lines');
h = 8000; % height (m)
l = 45000; % length
s = struct('boundary', [], 'lines', [], 'wells', []);
s.boundary = [0, 0;
              11493, 0;
              11500.86, 0; 
              l, 0;
              45000, 687.49;
              45000, 1769.31;
              45000, 2400;
              45000, 3350;
              45000, 3600;
              45000, 4500;
              l, h; ...
              0, h; ...
              0, 4014.68;
              0, 3212.54;
              0, 2988.4;
              0, 1856.99;
              0, 1352.38;
              0, 542.3;
              0, 0];
                          
% Wells (to have cells with correct centroids)
s.wells{1} = [13200, 2062];         % I1 (HW)
s.wells{2} = [12776, 2070];         % I2 (FW)
s.wells{3} = [20318, 2375];         % I3 (mid)

% Main fault lines
s.lines{1} = load(fullfile(pth, 'fault_CTl_coords.txt'));
s.lines{2} = load(fullfile(pth, 'fault_CT_coords.txt'));
s.lines{3} = load(fullfile(pth, 'fault_sinthetic_coords.txt'));

% Fault divisions (main units)
s.lines{4} = [11731.2 547.707;
              11742.47 547.85];
s.lines{5} = [12260 1354.98;
              12277.94 1352.38];
s.lines{6} = [12343.8 1451.3;
              12370.68 1450.28];
s.lines{7} = [12775.88 1871.46;
              12805.6 1871.67];
s.lines{8} = [12890.4 1978.9;
              12920.24 1973.54];
s.lines{9} = [14471.44 2990.1;
              14586.57 2988.4];
%s.lines{10} = [14985.197161 3214.629194;
%               15112.9      3212.54];

% Fault divisions (Amp b)
s.lines{10} = [12343.8 1451.3;
               12370.68 1450.28];
s.lines{11} = [12446.64 1551.3;
               12466.8 1550];
s.lines{12} = [12549.48 1651.3;
               12570    1650];
s.lines{13} = [12652.32 1751.3;
               12676.28 1750];
%s.lines{14} = [12755.1552785797 1851.303729;
%               12782.5617754334 1850];

% LM1
s.lines{14} = load(fullfile(pth, 'LM1_bot_left.txt'));
s.lines{15} = load(fullfile(pth, 'LM1_bot_right.txt'));

% Marg A
s.lines{16} = load(fullfile(pth, 'MargA_bot_left.txt'));
s.lines{17} = load(fullfile(pth, 'MargA_bot_mid.txt'));
s.lines{18} = load(fullfile(pth, 'MargA_bot_right.txt'));
s.lines{19} = load(fullfile(pth, 'MargA_top_left.txt'));
s.lines{20} = load(fullfile(pth, 'MargA_top_mid.txt'));
s.lines{21} = load(fullfile(pth, 'MargA_top_right.txt'));

% Amp B
s.lines{22} = load(fullfile(pth, 'Ampb_bot_left.txt'));
s.lines{23} = load(fullfile(pth, 'Ampb_bot_mid.txt'));
s.lines{24} = load(fullfile(pth, 'Ampb_bot_right.txt'));
s.lines{25} = load(fullfile(pth, 'Ampb_top_left.txt'));
s.lines{26} = load(fullfile(pth, 'Ampb_top_mid.txt'));
s.lines{27} = load(fullfile(pth, 'Ampb_top_right.txt'));

% Top Mio
s.lines{28} = load(fullfile(pth, 'Miocene_top_left.txt'));
s.lines{29} = load(fullfile(pth, 'Miocene_top_mid.txt'));
s.lines{30} = load(fullfile(pth, 'Miocene_top_right.txt'));

% Save data
save('gom_datapoints_noAmpDiv.mat', 's');
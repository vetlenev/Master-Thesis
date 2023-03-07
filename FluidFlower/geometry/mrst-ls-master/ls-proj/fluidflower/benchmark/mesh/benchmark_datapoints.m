%% Make datapoints for mesh conforming to boundary, fault, folded layers
%
% We draw the stratigraphy using paths in inkscape (model_lines.svg). 
% Then, we export nodes
% by selecting a given path and then doing extensions>export>export XY. Finally,
% we transform the y values using the actual y location (in mm) of the 
% leftmost node, since they are initially not absolute coordinates that we
% want. The extension export XY is external and was downloaded from here
% https://alpha.inkscape.org/vectors/www.inkscapeforum.com/viewtopicdb3f.html?t=8826
% and added to inkscape.
%

% Boundary
pth = fullfile(mrstPath('ls-proj'), 'fluidflower/benchmark/mesh/fromSvg');
topb = load(fullfile(pth, 'topBoundary.txt'));
h = 1340; % height (cm)
l = 2860; % length
to_m = 1e-3; % convert measuring/drawing units (cm) to m
stratiPoints = struct('boundary', [], 'lines', [], 'wells', []);
stratiPoints.boundary = to_m*[topb;
                              l, h; ...
                              0, h; ...
                              topb(1,:)];
                          
% Wells (to have cells with correct centroids)
stratiPoints.wells{1} = to_m*[930, h-330];          %I1
stratiPoints.wells{2} = to_m*[1730, h-730];         %I2
stratiPoints.wells{3} = to_m*[1530, h-530];         %P1
stratiPoints.wells{4} = to_m*[1730, h-1130];        %P2

% Lines (1 fault, 2 lines are shear zone boundaries)
fl_left = load(fullfile(pth, 'fault_top_left_left.txt'));
fl_right = load(fullfile(pth, 'fault_top_left_right.txt'));
fr_left = load(fullfile(pth, 'fault_top_right_left.txt'));
fr_right = load(fullfile(pth, 'fault_top_right_right.txt'));
fbl = load(fullfile(pth, 'fault_bot_left.txt'));
fbr = load(fullfile(pth, 'fault_bot_right.txt'));
fbm1 = load(fullfile(pth, 'fault_bot_midseg1.txt'));
fbm2 = load(fullfile(pth, 'fault_bot_midseg2.txt'));

% top left fault
stratiPoints.lines{1} = to_m*(fl_left);
stratiPoints.lines{2} = to_m*(fl_right);
stratiPoints.lines{3} = to_m*([fl_left(1,:); fl_right(end,:)]);
stratiPoints.lines{4} = to_m*([fl_left(end,:); fl_right(1,:)]);
stratiPoints.lines{5} = to_m*([1101.1428	524.80204
                               1083.7455	518.65724]); 
stratiPoints.lines{6}  = to_m*([1104.5987 541.49352
                               1089.7245 533.53736]);
                           
% top right fault
stratiPoints.lines{7} = to_m*(fr_left);
stratiPoints.lines{8} = to_m*(fr_right);
stratiPoints.lines{9}  = to_m*([fr_left(1,:); fr_right(end,:)]);
stratiPoints.lines{10}  = to_m*([fr_left(end,:); fr_right(1,:)]);


% bottom fault
stratiPoints.lines{11} = to_m*(fbl);
stratiPoints.lines{12} = to_m*(fbr);
stratiPoints.lines{13} = to_m*(fbm1);
stratiPoints.lines{14} = to_m*(fbm2);
stratiPoints.lines{15} = to_m*([408.76505	916.00982
                                462.21133	911.06111]);
stratiPoints.lines{16} = to_m*([408.1395	1042.0134
                                419.24694	1054.5911]);
stratiPoints.lines{17} = to_m*([fbl(1,:); fbr(1,:)]);

% Top left fault to middle contour
stratiPoints.lines{18} = to_m*([1100.8658	583.40473
                                1103.0688	576.68157]);
stratiPoints.lines{19} = to_m*([1112.7583	572.6018
                                1112.7393	583.31033]);

% Top left
tbl_top = load(fullfile(pth, 'tbl_top.txt'));
tbl_midtop = load(fullfile(pth, 'tbl_midtop.txt'));
tbl_midbot = load(fullfile(pth, 'tbl_midbot.txt'));
tbl_bot = load(fullfile(pth, 'tbl_bot.txt'));

stratiPoints.lines{20} = to_m*(tbl_top);
stratiPoints.lines{21} = to_m*(tbl_midtop);
stratiPoints.lines{22} = to_m*(tbl_midbot);
stratiPoints.lines{23} = to_m*(tbl_bot);

% Top middle
tbm_midtop = load(fullfile(pth, 'tbm_midtop.txt'));
tbm_midbot = load(fullfile(pth, 'tbm_midbot.txt'));
tbm_bot = load(fullfile(pth, 'tbm_bot.txt'));

stratiPoints.lines{24} = to_m*(tbm_midtop);
stratiPoints.lines{25} = to_m*(tbm_midbot);
stratiPoints.lines{26} = to_m*(tbm_bot);

% Top right
tbr_top = load(fullfile(pth, 'tbr_top.txt'));
tbr_midtop = load(fullfile(pth, 'tbr_midtop.txt'));
tbr_midbot = load(fullfile(pth, 'tbr_midbot.txt'));
tbr_bot = load(fullfile(pth, 'tbr_bot.txt'));
lens = load(fullfile(pth, 'lens_right.txt'));

stratiPoints.lines{27} = to_m*(tbr_top);
stratiPoints.lines{28} = to_m*(tbr_midtop);
stratiPoints.lines{29} = to_m*(tbr_midbot);
stratiPoints.lines{30} = to_m*(tbr_bot);
stratiPoints.lines{31} = to_m*(lens);

% Bot left
bl_top = load(fullfile(pth, 'bl_top.txt'));
bl_midtop = load(fullfile(pth, 'bl_midtop.txt'));
bl_midmid = load(fullfile(pth, 'bl_midmid.txt'));
bl_topseal = load(fullfile(pth, 'bl_topseal.txt'));
bl_botseal = load(fullfile(pth, 'bl_botseal.txt'));
bl_bot = load(fullfile(pth, 'bl_bot.txt'));

stratiPoints.lines{32} = to_m*(bl_top);
stratiPoints.lines{33} = to_m*(bl_midtop);
stratiPoints.lines{34} = to_m*(bl_midmid);
stratiPoints.lines{35} = to_m*(bl_topseal);
stratiPoints.lines{36} = to_m*(bl_botseal);
stratiPoints.lines{37} = to_m*(bl_bot);

% Bot mid
bm_top = load(fullfile(pth, 'bm_top.txt'));
bm_midtop = load(fullfile(pth, 'bm_midtop.txt'));

stratiPoints.lines{38} = to_m*(bm_top);
stratiPoints.lines{39} = to_m*(bm_midtop);

% Bot right
br_top = load(fullfile(pth, 'br_top.txt'));
br_topseal = load(fullfile(pth, 'br_topseal.txt'));
br_botseal = load(fullfile(pth, 'br_botseal.txt'));
br_bot = load(fullfile(pth, 'br_bot.txt'));

stratiPoints.lines{40} = to_m*(br_top);
stratiPoints.lines{41} = to_m*(br_topseal);
stratiPoints.lines{42} = to_m*(br_botseal);
stratiPoints.lines{43} = to_m*(br_bot);
                            
% Corrections
%stratiPoints.lines{15}(end-9,:) = stratiPoints.lines{1}(end,:);
%stratiPoints.lines{15}(end-8,:) = stratiPoints.lines{2}(end,:);
%stratiPoints.lines{15}(end-11,:) = stratiPoints.lines{9}(1,:);

% Save data
save('benchmark_datapoints_v2.mat', 'stratiPoints');
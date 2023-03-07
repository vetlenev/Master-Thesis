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

%dir = '/Users/lluis/Documents/MATLAB/mrst-dev/mrst-ls/ls-proj/fluidflower/medium/mesh/Bilbo/';
dir = 'C:/Users/lsalo/matlab/mrst-dev/mrst-ls/ls-proj/fluidflower/medium/mesh/Bilbo/';

% Boundary
zmin = 0.169;
h = 53.5-zmin; % height (cm)
l = 93.4;      % length
to_m = 1e-2;   % convert measuring/drawing units (cm) to m
stratiPoints = struct('boundary', [], 'lines', [], 'wells', []);
stratiPoints.boundary = load([dir 'boundary.txt']);
stratiPoints.boundary(:,2) = stratiPoints.boundary(:,2) - zmin;
stratiPoints.boundary = stratiPoints.boundary*to_m;

% transforms
yn = @(y) h - y;
                          
% Wells (to have cells with correct centroids)
stratiPoints.wells{1} = to_m*[76.7, yn(6.7)];
stratiPoints.wells{2} = to_m*[24.2, yn(24.2)];

% Bottom to top, left to right lines
stratiPoints.lines{1} = load([dir 'E1_HW.txt']);
stratiPoints.lines{2} = load([dir 'E1_FW.txt']);
stratiPoints.lines{3} = load([dir 'CbotESF.txt']);
stratiPoints.lines{4} = load([dir 'CtopESF.txt']);
stratiPoints.lines{5} = load([dir 'ESFmid_top.txt']);
stratiPoints.lines{6} = load([dir 'Cmid_HW.txt']);
stratiPoints.lines{7} = load([dir 'Cmid_FW.txt']);
stratiPoints.lines{8} = load([dir 'Ctop_HW.txt']);
stratiPoints.lines{9} = load([dir 'Ctop_FW.txt']);
stratiPoints.lines{10} = load([dir 'ESFtop_bot.txt']);
stratiPoints.lines{11} = [30.561121	26.485064
                          30.204693	24.457754];

for n=1:numel(stratiPoints.lines)
    stratiPoints.lines{n}(:,2) = stratiPoints.lines{n}(:,2) - zmin; 
    stratiPoints.lines{n} = stratiPoints.lines{n}*to_m;
end

% Save data
save('bilbo_datapoints.mat', 'stratiPoints');
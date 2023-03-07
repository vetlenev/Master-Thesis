function [plts, opt] = getSimOpts(inj_type, mesh_size, model_case, inj_mult, D)
%
%
%

%% Plots
plts.fig3D = @() figure('Position', [0, 0, 1300, 650]);
plts.setAxProps = @(ax) set(ax, 'View'              , [65, 20]         , ...
                           'PlotBoxAspectRatio', [4.40, 1.86, 1.00], ...
                           'Projection'        , 'Perspective'     , ...
                           'Box'               , 'on');
                       
%% Options
% Mesh
if inj_type == 1
    opt.mesh = {mesh_size, ['G_composite_cellSize' num2str(mesh_size) '.mat']};
elseif inj_type == 2
%     opt.mesh = {mesh_size, ['G_composite_cellSize' num2str(mesh_size) '.mat']};
    opt.mesh = {mesh_size, ['AC07G_composite_cellSize' num2str(mesh_size) '.mat']};
elseif inj_type == 3
    opt.mesh = {mesh_size, ['bilboG_composite_cellSize' num2str(mesh_size) '.mat']};
end

%BCs
if inj_type < 3
    opt.topboundary = 'straight'; % 'straight' or 'curved'
elseif inj_type == 3
    opt.topboundary = 'curved';
end
opt.bctype = 'cp';           % cp: constant pressure at top; pvm: pore volume mult

% Injection protocol
if inj_type == 1    % Albus, cycle 02 (the one used for matching)
    opt.inj_loc = [0.45, 0.47-.02; 0.897-.165, 0.47-.16]; 
    opt.wellno = 2;                                  % n injectors                    
    opt.rate = {[0.1 0.2 0.5 1 1.5 2]*inj_mult(1), 1:5;     % {mL/min, times (min);  inj 1
                [1.5 1 0.5 0.2 0.1 0]*inj_mult(1), 1:5;     %   "   ,   "   ; rampdown 1
                [0.1 0.2 0.5 1 1.5 2]*inj_mult(2), 1:5;     %  mL/min, times (min)}  inj 2
                [1.5 1 0.5 0.2 0.1 0]*inj_mult(2), ...      % " rampdonwn 2
                [34, 75, 60, 60, 60]};               %  seconds, rampdown 2
    opt.schedule = {0, 55, [69 11], 154, 48*60};     % [startI1, endI1, startI2, endI2, end sim]                                   
elseif inj_type == 2    % Albus, cycle 07 (only 1 injector at the bottom)
    opt.inj_loc = [0.45, 0.47-.02]; 
    opt.wellno = 2;                                  % n injectors                    
    opt.rate = {[0.1 0.2 0.5 1 1.5 2]*inj_mult(1), 1:5;     % {mL/min, times (min);  inj 1
                [1.5 1 0.5 0.2 0.1 0]*inj_mult(1), ...
                [28, 81, 60, 60, 60]};               
    opt.schedule = {0, [288 33], 24*60}; 
elseif inj_type == 3    % Bilbo, cycle 01 or 02
    opt.inj_loc = [0.7670, 0.4663; 0.2420 0.2913];
    opt.wellno  = 2;
    opt.rate = {[0.5 1 1.5 2]*inj_mult(1), 1:3;     % {mL/min, times (min);  inj 1
                [1.5 1 0.5 0]*inj_mult(1), 1:3;     %   "   ,   "   ; rampdown 1
                [0.5 1 1.5 2]*inj_mult(2), 1:3;     %  mL/min, times (min)}  inj 2
                [1.5 1 0.5 0]*inj_mult(2), 1:3};  % " rampdonwn 2
    opt.schedule = {0, [69 13], [72 5], [151 13], 24*60};
end

% Other model options
if model_case == 1
    opt.hyster = 'on';
    opt.vapoil = false;         % must be set, and must be true or false
    opt.D = D;               % diffusion coeff. [m2/s]
elseif model_case == 2
    opt.hyster = 'on';
    opt.vapoil = false;         % must be set, and must be true or false
    opt.D = D;    
elseif model_case == 3
    opt.hyster = 'on';
    opt.vapoil = false;         % must be set, and must be true or false
    opt.D = D; 
elseif model_case == 4
    
end

end

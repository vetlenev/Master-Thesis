function [state0, models, schedule, isFine, setZeroTrans] = setupSlopedGrid(faceConstraint, dims, sizes, varargin)
%SETUPHORIZONTALGRID Set up cartesian grid with horizontal impermeable
%layer in middle of domain.
%   Impermeable layer either to be represented as face constraint or as
%   cell constraint.
   
if nargin == 2
    sizes = [10*nx, 1*ny, 10*nz];       
elseif nargin < 2
    dims = [200, 1, 50];
    sizes = [800, 1, 250];
    if nargin == 0
        faceConstraint = true;
    end
end

if isempty(varargin)
    trans_multiplier = 1;
else
    trans_multiplier = varargin{1};
end

G = cartGrid(dims, sizes);
G0 = computeGeometry(G); % original grid before adding sealing layers/faces

nx = G.cartDims(1); ny = G.cartDims(2); nz = G.cartDims(3);

m = tand(4); % aquifer slope


%m = (k_range(2) - k_range(1)) / (i_range(2) - i_range(1));
% G.nodes.coords(:,3) = G.nodes.coords(:,3) + m*G.nodes.coords(:,1) ...
%                         + (nz/10)*sin(4*pi*G.nodes.coords(:,1)/max(G.nodes.coords(:,1))) ...
%                         + (exp(-3*G.nodes.coords(:,1)/max(G.nodes.coords(:,1))) * ...
%                                     (nz/2)) .* cos(pi/2+6*pi*G.nodes.coords(:,1)/max(G.nodes.coords(:,1))); % linear slope of aquifer

G = computeGeometry(G);

[ii, jj, kk] = gridLogicalIndices(G);

% perm = logNormLayers(G.cartDims);
% perm = perm./max(perm);
% perm = perm.*200*milli*darcy;
perm = repmat(100*milli*darcy, G.cells.num, 1);

setZeroTrans = zeros(G.faces.num, 1);
sealingCellsPerm = 0.1*milli*darcy;

%k_pos = [nz/6, nz/4, nz/2, 2*nz/3];
k_pos = [nz/4];
% i_range = [0, nx/2; ...
%            3*nx/4, nx; ...
%            nx/3, 2*nx/3; ...
%            0, nx];  
i_range = [3*nx/4, nx]; %...
            %3*nx/4, nx];
       
all_sealing_faces = [];
%all_sealing_cells = zeros(G.cells.num, 1);
all_sealing_cells = {};
       
if faceConstraint
    %k_range = [fix(nz/2), fix(nz/2)+1];
    for i=1:numel(k_pos)
        k_range = [fix(k_pos(i)), fix(k_pos(i))+i]; % +1
        [sealing_faces, sealingCells] = addConfiningLayers(G0, 'type', 'faces', 'i_range', [i_range(i,1), i_range(i,2)], 'k_range', k_range);                    
        all_sealing_faces = [all_sealing_faces; sealing_faces];
        all_sealing_cells = all_sealing_cells | sealingCells;       
    end
   
else % cell constraints
    %k_range = [fix(nz/2)-1, fix(nz/2)+2];
    for i=1:numel(k_pos)
        k_range = [fix(k_pos(i))-1, fix(k_pos(i))+2];
        [sealing_faces, sealingCells] = addConfiningLayers(G0, 'type', 'cells', 'i_range', [i_range(i,1), i_range(i,2)], 'k_range', k_range);       
        
        all_sealing_faces = [all_sealing_faces; sealing_faces];
        %all_sealing_cells = all_sealing_cells | sealingCells;
        all_sealing_cells = cat(2, all_sealing_cells, sealingCells);
        
        % assign low permeability to sealing layer
        nxi = numel(unique(ii(sealingCells)));
        nzi = numel(unique(kk(sealingCells)));
        perm(sealingCells) = logNormLayers([nxi, 1, nzi], repmat(sealingCellsPerm, nzi, 1)); % mean permeability of 0.1 mD 
        %perm(sealingCells) = sealingCellsPerm;
    end       
    
end

setZeroTrans(all_sealing_faces) = 1; 
setZeroTrans = logical(setZeroTrans); 

G = computeGeometry(G);
[x, z] = deal(G.cells.centroids(:,1), G.cells.centroids(:,3));

rock = makeRock(G, perm, 0.3);
% add compressibility ??
fluid = initSimpleADIFluid('phases', 'WG', ...
                           'mu', [8e-4 3e-5]*Pascal*second,...
                           'rho', [1100 700].* kilogram/meter^3, ... % simulate supercritical CO2
                           'pRef', 100*barsa);

tot_time = 500*year;
inj_stop = 0.1;
pv_rate = 0.05;

pv = poreVolume(G, rock);
inj_rate = pv_rate*sum(pv)/(inj_stop*tot_time); % inject pv_rate of total pore volume over injection time

% Specify well information
bc = [];
nearWell = false(G.cells.num, 1);

W = verticalWell([], G, rock, nx, 1, nz, ...
    'type', 'rate', ...  % inject at constant rate
    'val', inj_rate, ... % volumetric injection rate
    'comp_i', [0 1]);    % inject CO2, not water

bf = boundaryFaces(G);
bfx = G.faces.centroids(bf, 1);
bfz = G.faces.centroids(bf, 3);
lz = max(bfz(bfx == min(bfx)));
top_left_faces = bf(bfx == min(bfx) & bfz <= k_range(1)*lz/nz);
% 
% bc = addBC(bc, top_left_faces, 'pressure', ...
%             fluid.rhoWS.*G.faces.centroids(top_left_faces, 3)*norm(gravity), ...
%             'sat', [1,0]);
bc = pside(bc, G, 'West', 0*barsa, 'sat', [1, 0]);
bc.value = bc.value + fluid.rhoWS.*G.faces.centroids(bc.face, 3)*norm(gravity);

horzWellDistRate = 0.1; % ratio of total horizontal length considered "close" to well
vertWellDistRate = 0.1;
for i = 1:numel(W)
    c = W(i).cells(1);
    hdist = abs(x - x(c));
    vdist = abs(z - z(c));
    % store cells close to well, if fine-scale needed
    nearWell(hdist < fix(horzWellDistRate*max(x)) & ...
              vdist < fix(vertWellDistRate*max(z))) = true;
end
if ~isempty(bc)
    nearWell(sum(G.faces.neighbors(bc.face, :), 2)) = true;
end
   
nsteps_after_inj = 300;
dt = rampupTimesteps(inj_stop*tot_time, inj_stop*tot_time/200, 10);
dt = [dt; repmat((1-inj_stop)*tot_time/nsteps_after_inj, nsteps_after_inj, 1)];

times = cumsum(dt)/year();
n_steps = numel(dt);
[~, inj_stop] = min(abs(times - inj_stop*times(end)));

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);
schedule.control(2) = schedule.control(1); % new control for well shutoff
schedule.control(2).W.status = 0;
schedule.step.control(inj_stop:n_steps) = 2; % swap active well from first to second at halfway  

T = getFaceTransmissibility(G, rock);
% ---------
%T(setZeroTrans) = 0;
% ---------

isFine = struct;
isFine.sealing = all_sealing_cells; % no sealing cells if face constraint
isFine.well = nearWell;

model = TwoPhaseWaterGasModel(G, rock, fluid, 1, 1, 'useCNVConvergence', true);

model.operators.T = T(model.operators.internalConn); % update internal transmissibility
model.operators.T_all(model.operators.internalConn) = model.operators.T;

model_fine = model;

% Apply transmissibility multiplier to sealing faces
if ~isempty(all_sealing_faces)
    transMult = false(G.faces.num, 1);
    
    map = false(G.faces.num, 1);
    map(all_sealing_faces) = true; 
        
    transMult(map) = true;
    T(transMult) = T(transMult).*trans_multiplier;
    %model.operators.T = model.operators.T_all(model.operators.internalConn); 
end
    
model_fine.operators.T = T(model_fine.operators.internalConn); % update internal transmissibility
model_fine.operators.T_all(model_fine.operators.internalConn) = model_fine.operators.T; % update internal trans for full transmissibility matrix (boundary faces remains unchanged)

model.extraStateOutput = true;
model_fine.extraStateOutput = true;

models = struct;
models.original = model;
models.fine = model_fine;

state0 = initResSol(G, fluid.rhoWS*norm(gravity)*z, [1,0]);

end

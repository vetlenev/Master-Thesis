function [state0, models, schedule, isFine, setZeroTrans] = setupHorizontalGrid(faceConstraint, dims, sizes, varargin)
%SETUPHORIZONTALGRID Set up cartesian grid with horizontal impermeable
%layer in middle of domain.
%   Impermeable layer either to be represented as face constraint or as
%   cell constraint.
   
if nargin == 2
    nx = dims(1); ny = dims(2); nz = dims(3);
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
G = computeGeometry(G);

nx = G.cartDims(1); ny = G.cartDims(2); nz = G.cartDims(3);
gridSizes = max(G.faces.centroids);
lx = gridSizes(1); ly = gridSizes(2); lz = gridSizes(3);
[ii, jj, kk] = gridLogicalIndices(G);

% perm = logNormLayers(G.cartDims);
% perm = perm./max(perm);
% perm = perm.*200*milli*darcy;
perm = repmat(100*milli*darcy, G.cells.num, 1);

x_range = [0, lx/2]; % cover entire domain
setZeroTrans = zeros(G.faces.num, 1);
sealingCellsPerm = 0.1*milli*darcy;

if faceConstraint    
    k_range = [fix(nz/2), fix(nz/2+1)]; % single cell  
    
    sealingCells = zeros(G.cells.num, 1); 
    sealing_faces = addConfiningFaces(G, 'x_range', x_range, 'k_range', k_range);           
   
else % cell constraints      
    z_rate = 0.02;
    z_range = [fix(lz/2-lz*z_rate), fix(lz/2+lz*z_rate)]; % very small thickness
    
    [sealing_faces, sealingCells] = addConfiningCells(G, 'type', 'horizontal', 'x_range', x_range, 'z_range', z_range);       
    
    % assign low permeability to sealing layer
    nxi = numel(unique(ii(sealingCells)));
    nzi = numel(unique(kk(sealingCells)));
    perm(sealingCells) = logNormCells([nxi, 1, nzi], repmat(sealingCellsPerm, nzi, 1)); % mean permeability of 0.1 mD 
    
end

setZeroTrans(sealing_faces) = 1; 
setZeroTrans = logical(setZeroTrans); 

G = computeGeometry(G);
[x, z] = deal(G.cells.centroids(:,1), G.cells.centroids(:,3));

rock = makeRock(G, perm, 0.3);
% add compressibility ??
fluid = initSimpleADIFluid('phases', 'WG', ...
                           'mu', [8e-4 3e-5]*Pascal*second,...
                           'rho', [1100 700].* kilogram/meter^3, ... % simulate supercritical CO2
                           'n', [2,2], ...
                           'smin', [0, 0], ...
                           'pRef', 100*barsa);

tot_time = 1000*year;
inj_stop = 0.05;
pv_rate = 0.2;

pv = poreVolume(G, rock);
inj_rate = pv_rate*sum(pv)/(inj_stop*tot_time); % inject pv_rate of total pore volume over injection time

% Specify well information
bc = [];
nearWell = false(G.cells.num, 1);

W = verticalWell([], G, rock, 1, 1, nz, ...
    'type', 'rate', ...  % inject at constant rate
    'val', inj_rate, ... % volumetric injection rate
    'comp_i', [0 1]);    % inject CO2, not water

%bc = pside(bc, G, 'East', fluid.rhoWS*unique(z)*norm(gravity), 'sat', [1, 0]);
bf = boundaryFaces(G);
bfx = G.faces.centroids(bf, 1);
east_faces = bf(bfx == max(bfx));
bc = addBC(bc, east_faces, 'pressure', fluid.rhoWS*unique(z)*norm(gravity), 'sat', [1,0]);

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
   
nsteps_after_inj = 150;
dt = rampupTimesteps(inj_stop*tot_time, inj_stop*tot_time/200, 10);
dt_after_inj = rampupTimesteps((1-inj_stop)*tot_time, (1-inj_stop)*tot_time/nsteps_after_inj, 8);
%dt = [dt; repmat((1-inj_stop)*tot_time/nsteps_after_inj, nsteps_after_inj, 1)];
dt = [dt; dt_after_inj];

times = cumsum(dt)/year();
n_steps = numel(dt);
[~, inj_stop] = min(abs(times - inj_stop*times(end)));

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);
schedule.control(2) = schedule.control(1); % new control for well shutoff
schedule.control(2).W.status = 0;
schedule.step.control(inj_stop:n_steps) = 2; % swap active well from first to second at halfway  

T = getFaceTransmissibility(G, rock);
disp('min trans before')
disp(min(T))
disp('mean trans before')
disp(mean(T))

isFine = struct;
isFine.sealing = sealingCells; % no sealing cells if face constraint
isFine.well = nearWell;

model = TwoPhaseWaterGasModel(G, rock, fluid, 1, 1, 'useCNVConvergence', true);
model_fine = model;

% Apply transmissibility multiplier to sealing faces
if ~isempty(sealing_faces)
    transMult = false(G.faces.num, 1);
    
    map = false(G.faces.num, 1);
    map(sealing_faces) = true; 
        
    transMult(map) = true;
    T(transMult) = T(transMult).*trans_multiplier;
    %model.operators.T = model.operators.T_all(model.operators.internalConn); 
end

model_fine.operators.T = T(model_fine.operators.internalConn); % update internal transmissibility
model_fine.operators.T_all(model_fine.operators.internalConn) = model_fine.operators.T; % update internal trans for full transmissibility matrix (boundary faces remains unchanged)

model.extraStateOutput = true;
model_fine.extraStateOutput = true;

disp('min trans after')
disp(min(model.operators.T_all))
disp('mean trans after')
disp(mean(model.operators.T_all))

models = struct;
models.original = model;
models.fine = model_fine;

state0 = initResSol(G, fluid.rhoWS*norm(gravity)*z, [1,0]);

end

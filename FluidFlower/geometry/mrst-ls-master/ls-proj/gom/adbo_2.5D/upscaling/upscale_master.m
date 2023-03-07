%% Upscale k, Capillary Pressure and Relative Permeability
% This code takes a certain grid and rock structures (with poro and perm
% fields) plus "rock curves" for Pc and kr and upscales Pc and kr to the
% grid scale. 
% Multiple regions with different rock curves (e.g. sand and
% clay material) are included, and the option to use an externally-upscaled
% poro and perm (e.g. from PREDICT) when upscaling kr is given.
% Methods used are CL or invasion-percolation for Pc and VL and CL for kr. 
% See e.g. Li & Benson, AWR (2015) appendix B and C.
clear, close all force
mrstModule add ls-proj ls-utils mrst-gui mimetic upscaling incomp ...
               deckformat mpfa ad-props coarsegrid ad-blackoil ad-core ...
               linearsolvers
gravity reset off
%rng('default')

% Options
nsim = 1000;                    % fault realizations and upscaled curves
opt.nval = 32;                  % desired entries in output datafile
opt.fault =  'predict_3D';      % 'test' or 'predict'
opt.strati = 'GoM';             % only for predict
opt.window = 'famp6';           % for GoM: 'famp1' to 'famp6'
opt.clay_vcl = 0.7;             % vcl for clay-rich layers in caprock [0.4 to 0.7]
opt.theta = [30, 30];           % water-grains contact angle [sand, clay]
opt.dir = 'z';                  % 'x' (horizontal) or 'z' (vertical)
opt.pc_mode = 'inv-per';        % 'oper', 'inv-per', 'sim_inc', 'sim_bo'
opt.t = 1;                      % tau eq. (1) Yang et al., IJGGC 2013
opt.kr_mode = 'dyn';            % VL (viscous lim.), CL (capillary lim.), both, or 'dyn' for dynamic (high-rate)
makeplots = 1;                  %
opt = getUpsOpt(opt);           % get stratigraphy

pc = cell(nsim, 1);  
sg = cell(nsim, 1);  
rocks = cell(nsim, 1);
fluids = cell(nsim, 1);
faults = cell(nsim, 1);
sect = cell(1);
Us = cell(1);
tic
opt_fault = opt.fault;
parfor n=1:nsim
%for n=1
    % 1. Grid & Rock
    if strcmp(opt_fault, 'test')
        [G, rock, fault3D, ~, ~, ~, CG] = upscaleFaultPerm(opt, makeplots);  
    else
        [G, rock, fault3D, mySect, fault2D, Us_it, CG] = ...
                                          upscaleFaultPerm(opt, makeplots);
    end
    
    % 2. Fluids
    [rock, fluid, fault3D] = getFluidForUps(G, rock, fault3D, opt);
    
    % 3. Upscale Pc
    sg{n} = []; pc{n} = []; sgmax{n} = [];
    %[pc{n}, sg{n}, fluid, sgmax{n}] = upscalePcReg(G, fluid, rock, ...
    %                                            opt, makeplots);
   
    % 4. Upscale kr
    if ~strcmp(opt.kr_mode, 'dyn')
        [krw{n}, krn{n}] = upscaleKrReg(G, CG, rock, fluid, fault3D, ...
                                       sg{n}, pc{n}, sgmax{n}, opt, makeplots);
    end
    
    % Container variables
    Gs{n} = G;
    CGs{n} = CG;
    rocks{n} = rock;
    fluids{n} = fluid;
    faults{n} = fault3D;
    faultSections{n} = fault2D;
    sect{n} = mySect;
    if n == 1, Us{n} = Us_it; end
    
    % Loop status
    if mod(n, 10) == 0
       disp(['Simulation ' num2str(n) '/' num2str(nsim) ' completed.'])
    end
end
t1 = toc;

% 4. dynamic upscaling
% This needs to be done AFTER a full run of k to select the realization(s)
% that will be used for dynamic upscaling (liness 122-164)
if strcmp(opt.kr_mode, 'dyn')
    % For testing, find minimum ss for fastest simulation:
    %[~, id] =  min(cellfun(@(x) numel(x.SegLen), faults));
    %[~, id] =  max(cellfun(@(x) x.Perm(end), faults));
    
    % Get id from line 164
    vpar_fit = nan(numel(id), 3);
    mae_min = nan(numel(id), 1);
    sgof = nan(opt.nval*2+1, 4, numel(id));
    for n=1:numel(id)
        % Pc
        [pc, sg, fluid, sgmax] = upscalePcReg(Gs{id(n)}, fluids{id(n)}, ...
                                              rocks{id(n)}, opt, makeplots);
        % kr
        [krw, krn, sg2, vpar_fit(n,:), mae_min(n)] = upscaleKrReg(Gs{id(n)}, CGs{id(n)}, ...
                    rocks{id(n)}, fluid, faults{id(n)}, sg, pc, sgmax, ...
                    faultSections{id(n)}, sect{id(n)}, Us{1}, opt, makeplots);
                
        % Saturations compatible and save to SGOF-type table
        sg_fin = [0 1e-3 linspace(0.01, min(max(sg),max(sg2)), opt.nval-2)];
        pc_fin = interp1(sg, pc, sg_fin);
        krw_fin = interp1(1-sg2, krw, 1-sg_fin);
        krn_fin = interp1(sg2, krn, sg_fin);
        sgof(1:opt.nval,:,n) = [sg_fin', krn_fin', krw_fin', pc_fin'/barsa];
        
        % 5. Add "imbibition" curve (no hysteresis in fault)
        sgof(opt.nval+2:end,:,n) = [sg_fin', krn_fin', krw_fin', pc_fin'/barsa];
    end
end

% 6. Save data
fname = ['krPc_' opt.window '_k' opt.dir opt.dir '_theta' ...
         num2str(opt.theta(2)) '_PkzzBin_0.13'];
tab = table(sgof(:,1,1), sgof(:,2,1), sgof(:,3,1), sgof(:,4,1));
tab.Properties.VariableNames = {'SGAS' 'KRG' 'KRWG' 'PCOG'};
writetable(tab, [fname '.txt'], 'Delimiter', 'tab');

%% Visualize k distros and select indices for Pc, Kr
% Visualize stratigraphy and fault (thickness corresponding to 1st realization)
sect{1}.plotStrati(faults{1}.Thick, opt.fDip, 'm'); 

% Visualize fault materials
% Choice can be 'randm' (random), 'maxX', 'minX', 'maxZ' or 'minZ' or a
% specific perm value and direction (display the closest to that value)
plotId = selectSimId('randm', faults, nsim);                % simulation index
faults{plotId}.plotMaterials(faultSections{plotId}{1}, sect{plotId}, ...
                             'm', Us{1}) 

% Visualize upscaled permeability
plotUpscaledPerm(faults, 3, 'all')

% MANUALLY Select k values (after visualization above)
if opt.clay_vcl == 0.65
    if strcmp(opt.window, 'famp1')
        edg_log = [-1.67 -1.33 -0.67 -0.33;   % x_low1 x_high1 z_low1 z_high1 [mD]
                   -5.67 -5.33 -5 -4.67];     % x_low2 ...
    elseif strcmp(opt.window, 'famp2')
        edg_log = [-3.33 -3 -0.67 -0.33; 
                   -6 -5.67 -5.33 -5];      
    end
elseif opt.clay_vcl == 0.7
    if strcmp(opt.fault, 'predict')
        if strcmp(opt.window, 'famp1')
            edg_log = [-3.33 -3 -0.67 -0.33;   
                       -6 -5.67 -5.33 -5]; 
        elseif strcmp(opt.window, 'famp2')
            edg_log = [-3.33 -3 -0.67 -0.33;
                       -6 -5.67 -5.33 -5];
        end
    elseif strcmp(opt.fault, 'predict_3D')
        if strcmp(opt.window, 'famp1')
            edg_log = [-2.41 -2 -1.58 -1.16 -1.16 -0.75;    % x_low1 x_high1 y_low1 y_low1 z_low1 z_high1 [md]
                       -2.41 -2 -1.58 -1.16 -5.33 -4.91];   % x_low2 ..
        elseif strcmp(opt.window, 'famp2')
            edg_log = [-2.41 -2 -2 -1.58 -1.16 -0.75;
                       -6.16 -5.75 -5.33 -4.91 -5.33 -4.91];
        elseif strcmp(opt.window, 'famp3')
            edg_log = [-2 -1.58 -2 -1.58 -1.16 -0.75;
                       -6.16 -5.75 -5.33 -4.91 -5.33 -4.91];
        elseif strcmp(opt.window, 'famp4')
           edg_log = [-2 -1.58 -1.58 -1.16 -0.75 -0.33;
                       -6.16 -5.75 -5.33 -4.91 -5.33 -4.91];
        elseif strcmp(opt.window, 'famp5')
            edg_log = [-2 -1.58 -1.58 -1.16 -0.75 -0.33;
                       -5.75 -5.33 -5.33 -4.91 -5.33 -4.91];
        elseif strcmp(opt.window, 'famp6')
            edg_log = [-0.75 -0.33 0.08 0.5 0.5 0.92;
                       -5.75 -5.33 0.08 0.5 -4.91 -4.5];
        end
    end
end
         
% Find indices
ids = false(nsim, 2);
for n=1:nsim
    klog = log10(faults{n}.Perm/(milli*darcy));
    if klog(1) > edg_log(1,1) && klog(1) < edg_log(1,2) && ...
       klog(2) > edg_log(1,3) && klog(2) < edg_log(1,4) && ...
       klog(3) > edg_log(1,5) && klog(3) < edg_log(1,6)
        ids(n,1) = true;
    elseif klog(1) > edg_log(2,1) && klog(1) < edg_log(2,2) && ...
           klog(2) > edg_log(2,3) && klog(2) < edg_log(2,4) && ...
           klog(3) > edg_log(2,5) && klog(3) < edg_log(2,6)
        ids(n,2) = true;
    end
    %if all(~isnan(ids)), break, end
end
id_all = {find(ids(:,1)), find(ids(:,2))}; % for kr-pc upscaling
nSeg_id = cellfun(@(x) numel(x.SegLen), faults(id_all{1}));
[~, id_nSeg_min] = min(nSeg_id);
nSeg_id2 = cellfun(@(x) numel(x.SegLen), faults(id_all{2}));
[~, id_nSeg_min2] = min(nSeg_id2);
id = [id_all{1}(id_nSeg_min) id_all{2}(id_nSeg_min2)];

% Visualize upscaled Pc and Kr with selected values
plotPcKr(sg, pc, krw, krn, fluids, nsim, ids, opt.kr_mode)

close all force
S=rng;
save([opt.window '_k' opt.dir opt.dir '_theta' ...
      num2str(opt.theta(2)) '_' opt.pc_mode '_' opt.kr_mode '_vcl' ...
      num2str(opt.clay_vcl) '_N' num2str(nsim) '.mat'], '-v7.3')
%% Compute permeability of a fault with clay smears
%
% In the configuration used here, FW on the left and HW is on the right.
%
% The smears are consistent with (1) the number of sand and clay layers at
% each side of the fault; (2) the thickness of these layers; (3) fault 
% displacement and thickness; (4) the residual friction angle of the sands
% and clays; (5) the SSFc; and (6) the maximum length of a smear segment in
% each subdomain (to control whether smears are broken or not).
% Based on these inputs (1) to (6), the following are calculated: (I.) the
% probability of encountering a cell with smear in each smear subdomain;
% (II.) the length of each smear segment within each subdomain; and (III.)
% the thickness of each smear (and accordingly, the thickness of each
% subdomain).
%
%__________________________________________________________________________
clc, clear, close all
mrstModule add ls-proj ls-utils coarsegrid upscaling incomp mimetic mpfa ...
           streamlines mrst-gui
       
%% O. Organize variable values 
% Filename and output Figs
caseName       = 'nSimSmear_6bedsIrregAllClay_fT3fdip60_phiVar30_SSFcVar_LsVar_permVar';
toPlot         = 'permDistr';

% Input variables
changeClay     = 1;
faultDip       = 60;
%faultDip       = 60;
%faultThickness = [0.1 0.5 1 1.5 2 3 4 5 7 10];
faultThickness = 3;
%clayResidual    = 5:20;
clayResidual   = [10 5 8 12 3 16 10 5 8 4 3 5];
%sandResidual   = 25:35;
sandResidual   = 30;
%SSFc           = [1:0.25:4 5 7 10];
%SSFc           = 10;
SSFc           =[4 3 6 2 10 8 4 3 6 5 10 8];
%maxSmearSegLen = {2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 'Lsmear'};
maxSmearSegLen = {'Lsmear', 20, 'Lsmear', 10, 'Lsmear', 'Lsmear', 'Lsmear', 20, 'Lsmear', 10, 'Lsmear', 'Lsmear'};
%claySandPerms  = [20, 0.01, 0.1; 20, 1e-3 0.01; 20, 1e-4, 1e-3; 20, 1e-5, 1e-4]*milli*darcy;
%claySandPerms  = [20, 1e-4, 0.001]*milli*darcy;
sPerms = repmat([20, 20], 12, 1);          % 1x2 or N(clay layer) x 2
cPerms = [1e-4 1e-3; 1e-5 1e-4; 1e-2 1e-1; 1e-5 1e-4; 1e-4 1e-3; 1e-3 1e-2; 1e-4 1e-3; 1e-5 1e-4; 1e-2 1e-1; 1e-5 1e-4; 1e-4 1e-3; 1e-3 1e-2];      % 
[fdip, fT, phi_c, phi_s, SSFc, len, sandPerms, clayPerms, nSim] = smearParamVals(changeClay, ...
    faultDip, faultThickness, clayResidual, sandResidual, SSFc, maxSmearSegLen, sPerms, cPerms);
%sandPerms = [perms(:,1), perms(:,1)];   clayPerms = perms(:,2:3);
if changeClay == 1
   nDiffClays = size(SSFc, 2);
else
   nDiffClays = 1;
end


%% 1. Initialize
% Input independent variables
% Note: fw and hw layers' thickness (T) is the true layer thickness
%       to agree with Egholm et a. (2008). Relation to the apparent thickness
%       is Tap = T/(cosBl*sinB) where Bl is layer dip angle and B is
%       fault dip angle.
% Note: Collapse T and isclay in both fw and hw so that adjacent layers of
%       the same lithology are just one (eg [0 1] and not [0 1 1 1])

%nSimSmear = 100;
tic
%parfor n=1:nSim
for n=1:nSim
    %disp(['Simulation nr. = ' num2str(n) '/' num2str(nSim)]);
    f    = struct('T', fT(n), 'dip', fdip(n));     % fThickness not from grid
    sand = struct('perm', sandPerms((n-1)*nDiffClays+1:n*nDiffClays,:)*milli*darcy, 'poro', 0.3, 'phi', phi_s(n));
    clay = struct('perm', clayPerms((n-1)*nDiffClays+1:n*nDiffClays,:)*milli*darcy, 'poro', 0.1, 'phi', phi_c(n, :), 'SSFc', SSFc(n, :));
    fw   = struct('T', [10 30 20 10 20 10], 'isclay', repmat([true true], 1, 3), 'Bl', 0);
    hw   = struct('T', [5 20 25 25 10 15], 'isclay', repmat([true true], 1, 3), 'Bl', 0);
    g    = struct('resL', 1, 'resT', 0.1);          % grid resolution [m]
    smear.inL  = len(n, :);                         % Smear segment length [m] 
    smear.modeledAs = 'objects';                    % algorithm to place smears
    smear.tolerance = 0.01;                         % tolerance in P matching  
    displayFigures = [0 0 0];                       % [grid, smear location, perm] 

    %% 2. Create and populate grid
    % Cartesian grid in MRST has indexing that starts at (0, 0) coordinate
    % and runs faster on the X axis. Mapping grid and values for properties 
    % are stored in structure M.

    % 2.1 Compute initial values
    [f, sand, clay, fw, hw, smear] = smearOrganizeIn(f, sand, clay, fw, hw, smear);

    % 2.2 Create and plot grid
    nels = max([round(f.T/g.resT), round(f.D/g.resL)]);
    G = computeGeometry(cartGrid([nels, nels], [f.T f.D]));
    G.xzFaceDim = [f.T/nels, f.D/nels];

    % f1 = figure(1);
    % plotToolbar(G, G, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1); 
    % axis equal; colorbar; xlim([0 f.T]); ylim([0 f.D]); 
    % xlabel('x [m]'), ylabel('y [m]')
    % title(['Number of cells = ' num2str(G.cells.num)])
    % set(f1, 'position', [400, 0, 350, 900]); view([-30, 75]) 

    % 2.3 Mapping matrix structure
    M0 = smearMap(G, fw, hw, f, sand, clay, smear);

    % 2.4 Smears (final 1s)
    %     Continuity along possible diagonals 
    [M0, smear] = smearContinuity(fw, f, M0, sand, clay, smear);

    %     Smear distribution
    nSimSmear = [5 10 20 30 50 100 200 300 500 1000 2000 5000];
%     nSimSmear = [5 10 50];
%     P    = cell(numel(nSimSmear), 1);
%     poroMed = cell(numel(nSimSmear), 1);
%     permMed = cell(numel(nSimSmear), 1);
%     poroPctErr = cell(numel(nSimSmear), 1);
%     permPctErr = cell(numel(nSimSmear), 1);
    parfor k = 1:numel(nSimSmear)
%    for k = 1:numel(nSimSmear)    
        Pv    = cell(nSimSmear(k), 1);
        Porov = nan(nSimSmear(k),1);
        Permv = nan(nSimSmear(k),3);
        disp([num2str(k) '/' num2str(numel(nSimSmear))])
        valIds = 1:nSimSmear(k);
        for j = 1:nSimSmear(k)  % For single-parameter variation assessment, we need to avg smear distr. influence
            if any(M0.Psmear < 0.999) % place objects (segments) of Lsmear, and match P(smear)
                assert(strcmp(smear.modeledAs, 'objects'));
                verbose = 0;
                [M, Pv{j}] = smearObjects(M0, smear, f, G, verbose);
            elseif isempty(M0.Psmear)
                error('No smear: P(smear) = 0')
            else
                M = M0;
                assert(max(1 -M.Psmear) < 0.0011)
                Pv{j} = 1;
                disp('Continuous smears: P(smear) = 1')
            end            
            
            % assign porosity and permeability to each grid cell
            [poroG, permG] = poroPermToGrid(G, M, f, sand, clay);
            
            
            %% 3. Compute equivalent/upscaled permeability and porostiy
            % Porosity (additive)
            % Here, for 1 cell in upscaled grid and regular cells we could just do
            % Poro = (sum(sum(M.vals))*clay.poro +  ...
            %        (G.cartDims(1)^2 - sum(sum(M.vals)))*sand.poro) / G.cartDims(1)^2;
            p = partitionCartGrid(G.cartDims, [1 1]);
            CG = generateCoarseGrid(G, p);
            Porov(j) = accumarray(p, poroG)./accumarray(p,1);
            
            % Upscaled 2D Permeability tensor (Perm is not an additive property)
            gcellsCheck      = 0;               % check if Perm obtained with grid with AR of only 2 gives same Perm
            kbc              = 'sealed';        % 'sealed' or 'open' (periodic)
            kmethod          = 'mpfa';          % tpfa, mpfa or mimetic
            kacc             = 1;               % use mex acceleration in mpfa
            outfl            = 0;               % compare outflux of fine and coarse models
            outUps           = [1 10];          % [nx, ny] in coarse grid to compare outflow
            plotFigs         = 0;
            
            [K] = smearPermUpscaled(G, permG, f, 0, gcellsCheck, kbc, kmethod, kacc, outfl, outUps);
            %disp(['k_xx = ' num2str(K(1)/(milli*darcy)) '. k_zz = ' num2str(K(4)/(milli*darcy))])
            
            % Compute perm in third dimension and save the 3 independent components
            smearCells = sum(sum(M.vals));
            nct = G.cells.num;
            kyy = (smearCells/nct)*clay.perm(2) + ((nct - smearCells)/nct)*sand.perm(2);
            Permv(j,:) = [K(1), kyy, K(4)]./(milli*darcy);      % local fault perm [across, along-y, along-z]
            
            % Display progress
            if mod(j, 100) == 0
                disp([num2str(j) ' realizations out of ' num2str(nSimSmear) ...
                      '. Case ' num2str(n) ' out of ' num2str(nSim)]);
            end
            
            if max(1-M0.Psmear) < 0.0011
                valIds = 1;
                disp('M.Psmear > 1')
                break
            end
        end
        
        % Save poro and perm distributions (k loop)
        Poro(k).all              = Porov(valIds, 1);
        Perm(k).all              = Permv(valIds, :);
        
        %Pv = cell2mat(Pv);
        %bePoro              = median(Porov(valIds));
        %Poro(n).estimate    = bePoro;
        %Poro(n).pctErr      = [abs(bePoro - prctile(Porov(valIds), 10, 1)), ...
        %                       abs(bePoro - prctile(Porov(valIds), 90, 1))];
        
        
        %bePerm              = median(Permv(valIds, :), 1);
        %Perm(n).estimate    = bePerm;
        %Perm(n).pctErr      = [abs(bePerm - prctile(Permv(valIds, :), 10, 1)), ...
        %                       abs(bePerm - prctile(Permv(valIds, :), 90, 1))];
        %P                   = mean(Pv(valIds, :));

%         poroMed{k}(n)      = median(Porov(valIds));
%         poroPctErr{k}(n,:) = [abs(poroMed{k}(n) - prctile(Porov(valIds), 10, 1)), ...
%                               abs(poroMed{k}(n) - prctile(Porov(valIds), 90, 1))];
%         permMed{k}(n, :)   = median(Permv(valIds, :), 1);
%         permPctErr{k}(n,:) = [abs(permMed{k}(n, :) - prctile(Permv(valIds, :), 10, 1)), ...
%                               abs(permMed{k}(n, :) - prctile(Permv(valIds, :), 90, 1))];
%         P{k}               = mean(Pv(valIds, :));

    end
    
    % Save poro and perm distributions (n loop)
%     Poro(n).all              = Porov(valIds, 1);
%     Perm(n).all              = Permv(valIds, :);
end
toc


%% Single-parameter variation analysis
% if exist('poroMed', 'var') && iscell(poroMed)
%    Poro.estimate = poroMed;
%    Poro.pctErr = poroPctErr;
%    Perm.estimate = permMed;
%    Perm.pctErr = permPctErr;
% end

% Analysis
for n=1:numel(Perm)
    Perm(n).med = median(Perm(n).all, 1);
    Perm(n).pctErr = [abs(Perm(n).med - prctile(Perm(n).all, 10, 1)), ...
                      abs(Perm(n).med - prctile(Perm(n).all, 90, 1))];
end

% Plots
claySandPerms = [max(sPerms(:, 1)), min(cPerms(:, 1)), min(cPerms(:, 2))]*milli*darcy;
singleParameterPlots(toPlot, nSimSmear, claySandPerms, Perm, Poro, clayResidual, ...
                     sandResidual, SSFc, len, faultDip, faultThickness);
% 
% 
disp('Data will be saved. CHECK filename and press enter.')
pause
close all force
save([caseName '.mat'])
disp(['Data saved in ' pwd])


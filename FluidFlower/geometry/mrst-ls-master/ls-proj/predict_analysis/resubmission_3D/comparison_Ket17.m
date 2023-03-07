%% Comparison of PREDICT output with experiments in Kettermann et al., JGR:SE (2017)
% 
%
clear
close all force
%rng('default')          % check repeatability

%% 1. Load Required MRST Modules
mrstModule add mrst-gui coarsegrid upscaling incomp mpfa mimetic
mrstVerbose on     

%% 2. Define Model and Upscale Permeability
% 2.1 Mandatory input parameters
thickness = {[0.05 0.01], 0.06};                                     % [m]
vcl       = {[0 1], 0};                                                 % fraction [-]
%dip       = [0, 0];                                                        % [deg.]
faultDip  = 60;                                                             % [deg.]
zf        = [0, 0];                                                         % [FW, HW], [m]
zmax      = {repelem(1, numel(vcl{1})), repelem(1, numel(vcl{2}))};         % {FW, HW}
dim       = 3;                    % dimensions (2 = 2D, 3 = 3D)

% 2.2 Optional input parameters
%     We indicate the measured layer permeabilties and set failure type as
%     hybrid.
maxPerm = 10^-4*10^-3/(998.23*9.81)/(milli*darcy);      % k = K*mu/(rho*g) (sand) [mD]
clayPerm = 4*10^-8*10^-3/(998.23*9.81)/(milli*darcy);   % [mD]
FW_perm_vals = [nan, clayPerm];        
clayFailure = 'hybrid';                               % default (not indicated) is shear
rho = 0.6;                                            % Corr. coeff. for multivariate distributions

% 2.3 Flow upscaling options and number of simulations
U.useAcceleration = 1;          % 1 requires MEX setup, 0 otherwise (slower for MPFA).
U.method          = 'tpfa';     % 'tpfa' recommended for 3D
U.coarseDims      = [1 1 1];    % desired n cells [x, y, z] in coarse grid
U.flexible        = true;
Nsim              = 100;        % Number of 3D simulations/realizations

% 2.4 Define Stratigraphy and FaultedSection objects
footwall = Stratigraphy(thickness{1}, vcl{1}, ...
                        'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'IsHW', 1, ...
                           'NumLayersFW', footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', zmax{2});
footwall.Perm = FW_perm_vals;

% Instantiate FaultedSection object (Strati in Faulted Section)
mySect = FaultedSection(footwall, hangingwall, faultDip, ...
                        'maxPerm', maxPerm, 'failure', clayFailure);

% 2.5 Get material distributions
mySect = mySect.getMatPropDistr();
% Get along-strike segmentation
nSeg = getNSeg(mySect.Vcl, mySect.IsClayVcl, mySect.DepthFaulting);

% 2.6 Generate intermediate variable samples, calculate smear dimensions 
%     and upscale permeability.
assert(dim==3);
faultSections = cell(Nsim, 1);
smears = cell(Nsim, 1);
faults = cell(Nsim, 1);
upscaledPerm = zeros(Nsim, 3);
D = sum(mySect.Tap(mySect.FW.Id));      % displacement
nSeg_fcn = nSeg.fcn;
tstart = tic;
parfor n=1:Nsim    % parfor allowed if you have the parallel computing toolbox
%for n=1
    % Instantiate fault section and get segmentation for this realization
    myFaultSection = Fault2D(mySect, faultDip);
    myFault = Fault3D(myFaultSection, mySect);
    myFault = myFault.getSegmentationLength(U, nSeg_fcn);
    G = [];
    for k=1:numel(myFault.SegLen)
        % Get material property (intermediate variable) samples, and fix
        % fault thickness of current realization (3D only).
        myFaultSection = myFaultSection.getMaterialProperties(mySect, 'corrCoef', rho);
        myFaultSection.MatProps.thick = myFault.Thick;
        if isempty(G)
            G = makeFaultGrid(myFault.Thick, myFault.Disp, ...
                              myFault.Length, myFault.SegLen, U);
        end
        
        % Generate smear object with T, Tap, L, Lmax
        smear = Smear(mySect, myFaultSection, G, 1);
        
        % Place fault materials and assign cell-based properties in 2D section
        myFaultSection = myFaultSection.placeMaterials(mySect, smear, G);
        
        % Extrude 2D section to fill current segment
        myFault = myFault.assignExtrudedVals(G, myFaultSection, k);
        
        % Save results
        faultSections{n}{k} = myFaultSection;
        smears{n}{k} = smear;
    end

    % Compute 3D upscaled permeability distribution
    myFault = myFault.upscaleProps(G, U);
    
    % Save results
    faults{n} = myFault;
    if mod(n, 50) == 0
        disp(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
    end
end
telapsed = toc(tstart);

%% 3. Output Analysis
% % 3.1 Visualize stratigraphy and fault (with thickness corresponding to 1st realization)
mySect.plotStrati(faults{1}.Thick, faultDip, 'cm');  
% 
% % 3.2 Visualize intermediate variables
% % We define a given parent material (id from 1 to n of materials in stratigraphy), 
% % and generate histograms and correlation matrix plots.
% layerId = 2;                                            
% plotMatPropsHist(faultSections, smears, mySect, layerId, dim) 
% % MatProps correlations
% [R, P] = plotMatPropsCorr(faultSections, mySect, layerId, dim);
% if dim==3
%     plotSeg(faults, nSeg)
% end

% 3.3 Visualize fault materials
% Visualization for one realization. Choice can be 'randm' (random), 'maxX' 
% (realization with maximum upscaled permeability in across the fault), 'minX', 
% 'maxZ' or 'minZ'.
% General fault materials and perm view
plotId = selectSimId('randm', faults, Nsim);                % simulation index
%plotId = 1;
faults{plotId}.plotMaterials(faultSections{1}{1}, mySect, ...
                            'cm', U) 

% 3.4. Visualize upscaled permeability
% Plot upscaled permeability distributions (all simulations)
plotUpscaledPerm(faults, dim)


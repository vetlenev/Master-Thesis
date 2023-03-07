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
thick_all = {{[0.01 0.01], 0.02}, {[0.03 0.01], 0.04}, ...
             {[0.05 0.01], 0.06}, {[0.07 0.01], 0.08}, ...
             {[0.005 0.04 0.005], 0.05}, {[0.02 0.005 0.04 0.005], 0.07} ...
             {[0.04 0.005 0.04 0.005], 0.09}};         
vcl_all    = {{[0 1], 0}, {[1 0 1], 0}, {[0 1 0 1], 0}};                                                    
%dip       = [0, 0];                                                        % [deg.]
faultDip  = 60;                                                             % [deg.]
zf        = [0, 0];                                                         % [FW, HW], [m]
dim       = 3;                    % dimensions (2 = 2D, 3 = 3D)

% Optional input parameters
%     We indicate the measured perms and set failure type as hybrid.
maxPerm = 10^-4*10^-3/(998.23*9.81)/(milli*darcy);      % k = K*mu/(rho*g) (sand) [mD]
clayPerm = 4*10^-8*10^-3/(998.23*9.81)/(milli*darcy);   % [mD]
clayFailure = 'hybrid';                                 % default (not indicated) is shear
rho = 0.6;                                              % Corr. coeff. for multivariate distributions
         
% Flow upscaling options and number of simulations
U.useAcceleration = 1;           % 1 requires MEX setup, 0 otherwise (slower for MPFA).
U.method          = 'tpfa';      % 'tpfa' recommended for 3D
U.coarseDims      = [1 1 1];     % desired n cells [x, y, z] in coarse grid
Nsim              = 1000;        % Number of 3D simulations/realizations

SSF = [2 4 6 8 5 7 9];
faults = cell(Nsim, numel(SSF));
assert(dim==3);
tstart = tic;
for j=1:numel(SSF)
thickness = thick_all{j};                                                   % [m]
if j < 5
    vcl = vcl_all{1};
elseif j==5
    vcl = vcl_all{2};
else
    vcl = vcl_all{3};
end
zmax = {repelem(1, numel(vcl{1})), repelem(1, numel(vcl{2}))};              % {FW, HW}

% Define Stratigraphy and FaultedSection objects
footwall = Stratigraphy(thickness{1}, vcl{1}, ...
                        'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'IsHW', 1, ...
                           'NumLayersFW', footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', zmax{2});
if j < 5
    footwall.Perm = [nan, clayPerm];
elseif j==5
    footwall.Perm = [clayPerm, nan, clayPerm];
else
    footwall.Perm = [nan, clayPerm, nan, clayPerm];
end

% Instantiate FaultedSection object (Strati in Faulted Section)
mySect = FaultedSection(footwall, hangingwall, faultDip, ...
                        'maxPerm', maxPerm, 'failure', clayFailure);

% Get material distributions
mySect = mySect.getMatPropDistr();
% Get along-strike segmentation
nSeg = getNSeg(mySect.Vcl, mySect.IsClayVcl, mySect.DepthFaulting);

% Generate intermediate variable samples, calculate smear dimensions 
% and upscale permeability.
nSeg_fcn = nSeg.fcn;
parfor n=1:Nsim    
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
    end

    % Compute 3D upscaled permeability distribution
    myFault = myFault.upscaleProps(G, U);
    
    % Save results
    faults{n, j} = myFault;
    if mod(n, 50) == 0
        disp(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
    end
end

end
telapsed = toc(tstart);

% Comparison figure
comparisonPlot(faults, SSF)

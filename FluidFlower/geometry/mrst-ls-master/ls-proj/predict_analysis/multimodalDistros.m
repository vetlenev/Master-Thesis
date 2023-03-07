%% 1. Clay Smear Modeling Leads to Multimodal Fault Permeability Distributions
%
% In this script, we illustrate how stratigraphies with different clay
% contents lead to multimodal fault permeability distributions.
% 

clear
close all force

%% Mrst modules. 
% Run startup (enough for grid) + add required for flow upscaling. With 
% either development or release version of MRST.
%
% We use the following MRST utilities:
%   * merge_options.m -->
%   *
%   * 
mrstModule add mrst-gui coarsegrid upscaling incomp mpfa


%% Define model and upscale permeability
idSeq = [0.2, 0.4, 0.6, 0.8] ;            % strati identifiers, smear source fraction in seq.
vclSeq = [0.205 0.265 0.3575 0.5275];     % avg vcl sequence
fname = [];         % [] (empty) to not save data.
dim = 2;

% Mandatory Input parameters
thickness = {{[20 10 20 10 40], [5 20 5 20 5 20 5 20]}, ...
             {[10 20 10 20 10 10 10 10], [5 10 15 30 20 20]}, ...
             {[20 10 20 10 20 20], [15 20 30 20 15]}, ...
             {[10 10 20 20 20 20], [20 5 15 10 20 30]}};
vcl       = {{[0.2 0.4 0.1 0.5 0.1], [0.5 0.2 0.4 0.1 0.5 0.3 0.4 0.05]}, ...
             {[0.4 0.1 0.5 0.2 0.6 0.1 0.5 0.2], [0.6 0.05 0.5 0.1 0.4 0.1]}, ...
             {[0.6 0.1 0.5 0.2 0.4 0.05], [0.4 0.1 0.5 0.2 0.7]}, ...
             {[0.6 0.4 0.1 0.6 0.8 0.7], [0.2, 0.7, 0.8, 0.4, 0.5, 0.6]}};
faultDip  = 70;

% Optional Input parameters
nl   = [numel(vcl{1}{1}), numel(vcl{1}{2}); ...
        numel(vcl{2}{1}), numel(vcl{2}{2}); ...
        numel(vcl{3}{1}), numel(vcl{3}{2}); ...
        numel(vcl{4}{1}), numel(vcl{4}{2})];  % just for convenience here
zf   = [1000, 1000];    % m
zmax = {{repelem(2000, 1, nl(1,1)), repelem(2000, 1, nl(1,2))}, ...
        {repelem(2000, 1, nl(2,1)), repelem(2000, 1, nl(2,2))}, ...
        {repelem(2000, 1, nl(3,1)); repelem(2000, 1, nl(3,2))}, ...
        {repelem(2000, 1, nl(4,1)); repelem(2000, 1, nl(4,2))}};
cm = 'kao';                                     % predominant clay mineral
maxPerm = 1000;                                 % cap max perm? [mD]
rho = 0.6;                                      % Corr. coeff. for multiv. distr.    

% Flow upscaling options
U.useAcceleration = 1;          % requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm
U.coarseDims      = [1 1 1];

% Prepare loop
nStrat = numel(thickness);
Nsim = 1000;
faults_all = cell(nStrat, 1);
smears_all = cell(nStrat, 1);
sect_all = cell(nStrat, 1);
C = zeros(Nsim,nStrat,2);
Cx = nan(Nsim, 1); Cz = nan(Nsim, 1);
tic
for j=1:nStrat           % For each stratigraphy
    disp(['Stratigraphy ' num2str(j) ' / ' num2str(nStrat) ' in progress...'])
    
    % FW and HW
    footwall = Stratigraphy(thickness{j}{1}, vcl{j}{1}, ...
                            'DepthFaulting', zf(1), ...
                            'DepthBurial', zmax{j}{1}, 'ClayMine', cm);
    hangingwall = Stratigraphy(thickness{j}{2}, vcl{j}{2}, 'IsHW', 1, ...
                               'NumLayersFW', footwall.NumLayers, ...
                               'DepthFaulting', zf(2), ...
                               'DepthBurial', zmax{j}{2}, 'ClayMine', cm);
    
    % Strati in Faulted Section
    mySect = FaultedSection(footwall, hangingwall, faultDip, ...
                            'maxPerm', maxPerm);
    
    % Get material property distributions
    mySect = mySect.getMatPropDistr();
    
    % Base grid
    D = sum(mySect.Tap(mySect.FW.Id));
    T0 = 1;
    G0 = makeFaultGrid(T0, D);
    
    % Loop over realizations
    faults = cell(Nsim, 1);
    smears = cell(Nsim, 1);
    parfor n=1:Nsim  
        % Generate fault object with properties for each realization
        myFault = Fault(mySect, faultDip);
        
        % Get dependent variables
        myFault = myFault.getMaterialProperties(mySect, 'corrCoef', rho);
        
        % Update grid dimensions with sampled fault thickness
        G = updateGrid(G0, myFault.MatProps.thick);
        
        % Generate smear object with T, Tap, L, Lmax
        smear = Smear(mySect, myFault, G, 1);
        
        % Place fault materials and assign cell-based properties
        myFault = myFault.placeMaterials(mySect, smear, G);
        
        % Compute upscaled permeability distribution
        myFault = myFault.upscaleProps(G, U);
        
        % Find boundary to boundary connections
        c = findSandConn(myFault.MatMap.vals);
        
        % Save result
        faults{n} = myFault;
        smears{n} = smear;
        Cx(n) = c.bc(1); Cz(n) = c.bc(2);
        if mod(n, 500) == 0
            disp([num2str(n) ' realizations out of ' num2str(Nsim), ...
                ' completed.'])
        end
    end
    sect_all{j} = mySect;
    faults_all{j} = faults;
    smears_all{j} = smears;
    C(:,j,1) = Cx;  C(:,j,2) = Cz;
    
    disp(['Stratigraphy ' num2str(j) ' finished.'])
    disp('***************************************************************')
end
telapsed = toc;

% Save data?
if ~isempty(fname)
    disp(['ATTENTION: data saved in: ' pwd ' with filename ' fname])    
    save([fname '.mat'],'-v7.3') % larger than 2GB
end


%% Output analysis
stratId = 3;

% Visualize Strati, fault thickness of 1st realization
sect_all{stratId}.plotStrati(faults_all{stratId}{1}.MatProps.thick, faultDip, 'm'); 

% Histograms for each MatProp (all sims, we select one stratigraphic layer)
% This should plot for all realizations that contain the given id.
%layerId = 2;                                            
%plotMatPropsHist(faults_all{stratId}, smears_all{stratId}, ...
%                 sect_all{stratId}, layerId, dim) 

% MatProps correlations
%[R, P] = plotMatPropsCorr(faults_all{stratId}, sect_all{stratId}, 2, dim);

% General fault materials and perm view
plotId = selectSimId('randm', faults_all{stratId}, Nsim);                % simulation index
faults_all{stratId}{plotId}.plotMaterials(sect_all{stratId}, G0) 

% Plot upscaled Poro and Perm (all sims, 3 directions)
plotUpscaledPerm(faults_all{stratId}, dim, 'histOnly')

% Plot connected sand pathways
plotConnected(C, vclSeq)

%% 1. Clay Smear Modeling Leads to Multimodal Fault Permeability Distributions
%
% Multiple vcl and comparison with Spe02, Grant.
%
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
idSeq = {'0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', '0.45', ...
         '0.5', '0.55', '0.6', '0.65', '0.7', '0.75', '0.8', '0.85'};       % strati identifiers, smear source fraction in seq.
fname = 'clayContentAndFperm2_Nsim1000_8x2';    % [] (empty) to not save data.

[thickness, vcl, faultDip, zf, zmax, cm] = getInputs(fname);
maxPerm = [];                                   % cap max perm? [mD]
rho = 0.6;                                      % Corr. coeff. for multiv. distr.    

% Flow upscaling options
U.useAcceleration = 1;          % requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm

% Prepare loop
nStrat = numel(thickness);
Nsim = 1000;
faults_all = cell(nStrat, 1);
smears_all = cell(nStrat, 1);
sect_all = cell(nStrat, 1);
tic
for j=1:nStrat           % For each stratigraphy
    disp(['Stratigraphy ' num2str(j) ' / ' num2str(nStrat) ' in progress...'])
    
    % FW and HW
    footwall = Stratigraphy(thickness{j}{1}, vcl{j}{1}, ...
                            'DepthFaulting', zf(j,1), ...
                            'DepthBurial', zmax{j}{1}, 'ClayMine', cm(j, :));
    hangingwall = Stratigraphy(thickness{j}{2}, vcl{j}{2}, 'IsHW', 1, ...
                               'NumLayersFW', footwall.NumLayers, ...
                               'DepthFaulting', zf(j,2), ...
                               'DepthBurial', zmax{j}{2}, 'ClayMine', cm(j, :));
    
    % Strati in Faulted Section
    mySect = FaultedSection(footwall, hangingwall, faultDip, ...
                            'maxPerm', maxPerm);
    
    % Get material property distributions
    mySect = mySect.getMatPropDistr();
    
    % Loop over realizations
    faults = cell(Nsim, 1);
    smears = cell(Nsim, 1);
    parfor n=1:Nsim  
        % Generate fault object with properties for each realization
        myFault = Fault(mySect, faultDip);
        
        % Get dependent variables
        myFault = myFault.getMaterialProperties(mySect, 'corrCoef', rho);
        
        % Generate smear object with T, Tap, L, Lmax
        smear = Smear(mySect, myFault, 1);
        
        % Compute upscaled permeability distribution
        myFault = myFault.upscaleSmearPerm(mySect, smear, U);
        
        % Save result
        faults{n} = myFault;
        smears{n} = smear;        
        if mod(n, 500) == 0
            disp([num2str(n) ' realizations out of ' num2str(Nsim), ...
                ' completed.'])
        end
    end
    sect_all{j} = mySect;
    faults_all{j} = faults;
    smears_all{j} = smears;
    
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
stratId = 14;

% Visualize Strati, fault thickness of 1st realization
sect_all{stratId}.plotStrati(faults_all{stratId}{1}.MatProps.thick, faultDip, 'large'); 

% Histograms for each MatProp (all sims, we select one stratigraphic layer)
% This should plot for all realizations that contain the given id.
%layerId = 2;                                            
%plotMatPropsHist(faults_all{stratId}, smears_all{stratId}, ...
%                 sect_all{stratId}, layerId) 

% MatProps correlations
%[R, P] = plotMatPropsCorr(faults_all{stratId}, sect_all{stratId}, 2);

% General fault materials and perm view
%plotId = selectSimId('randm', faults_all{stratId}, Nsim);                % simulation index
%faults_all{stratId}{plotId}.plotMaterials(sect_all{stratId}) 

% Plot upscaled Poro and Perm (all sims, 3 directions)
%plotUpscaledPerm(faults_all{stratId}, 'histOnly')

% Boxplot graph
%plotMultipleBoxPlot(faults_all, Nsim, idSeq, thickness, vcl, ...
%                    zf, zmax, 'full')
plotHistComparison(faults_all{stratId}, thickness{stratId}, vcl{stratId}, ...
                   zf(stratId, :), zmax{stratId})
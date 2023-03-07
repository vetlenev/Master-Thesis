%% 1. Clay Smear Modeling Leads to Multimodal Fault Permeability Distributions
%
% In this script, we illustrate how 3 different stratigraphies with similar
% clay content lead to different permeability distributions, as a result of
% clay smear location and properties in the fault.
% The 3 stratigraphies are representative of the Miocene in the GoM, the
% Jurassic (Statfjord Fm.) in the northern North Sea, and Miocene (Mount 
% Messenger Fm.) in Taranaki basins.
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
idSeq = ['GoM'; 'NNS'; 'Tar'];           % strati identifiers, smear source fraction in seq.
fname = 'stratiControls_Nsim5k_3D';     % [] (empty) to not save data.

% Mandatory Input parameters
thickness = {{[6 3 1], [4 4 2]}, ...
             {[1 2 2 1.5 1 2.5], [1.5 2 2 1 1.5 2]}, ...
             {[0.5 1 .5 1 1 .5 .5 1 1 .5 .5 2], [.5 1 .5 1 1 .5 .5 2 1 .5 .5 1]}};
vcl       = {{[0.15 0.6 0.2], [0.15 0.5 0.25]}, ...
             {[0.5 0.1 0.4 0.2 0.5 0.3], [0.4 0.2 0.4 0.1 0.5 0.15]}, ...
             {[0.5 0.25 0.4 0.2 0.5 0.25 0.4 0.15 0.4 0.25 0.4 0.2], ...
              [0.4 0.225 0.4 0.2 0.4 0.25 0.4 0.2 0.5 0.2 0.5 0.2]}};
faultDip  = 70;
dim = 3;
unit_plot = 'm';

idl = 3;
vcl_sect = sum([thickness{idl}{1}.*vcl{idl}{1} thickness{idl}{2}.*vcl{idl}{2}]) / ...
           sum([thickness{idl}{1} thickness{idl}{2}]);   % check same section vcl

% Optional Input parameters
nl   = [numel(vcl{1}{1}), numel(vcl{1}{2}); ...
        numel(vcl{2}{1}), numel(vcl{2}{2}); ...
        numel(vcl{3}{1}), numel(vcl{3}{2})];  % just for convenience here
zf   = [200, 200; 700, 700; 1500, 1500];      % m
zmax = {{repelem(2000, 1, nl(1,1)), repelem(2000, 1, nl(1,2))}, ...
        {repelem(2000, 1, nl(2,1)), repelem(2000, 1, nl(2,2))}, ...
        {repelem(1500, 1, nl(3,1)); repelem(1500, 1, nl(3,2))}};
cm = 'kao';                                     % predominant clay mineral
maxPerm = 5000;                                 % cap max perm? [mD]
rho = 0.6;                                      % Corr. coeff. for multiv. distr.    

% Flow upscaling options
U.useAcceleration = 1;          % requires MEX and AMGCL setup
U.method          = 'tpfa';     % 'tpfa' in 3D
U.flexible        = false;
U.coarseDims      = [2 10 10];

% Prepare loop
nStrat = numel(thickness);
Nsim = 1000;
faults = cell(Nsim, nStrat);
faultSections_all = cell(1, Nsim);
%smears_all = cell(nStrat, 1);
sect_all = cell(nStrat, 1);
tic
for j=1:nStrat           % For each stratigraphy
    disp(['Stratigraphy ' num2str(j) ' / ' num2str(nStrat) ' in progress...'])
    
    % FW and HW
    footwall = Stratigraphy(thickness{j}{1}, vcl{j}{1}, ...
                            'DepthFaulting', zf(j, 1), ...
                            'DepthBurial', zmax{j}{1}, 'ClayMine', cm);
    hangingwall = Stratigraphy(thickness{j}{2}, vcl{j}{2}, 'IsHW', 1, ...
                               'NumLayersFW', footwall.NumLayers, ...
                               'DepthFaulting', zf(j, 2), ...
                               'DepthBurial', zmax{j}{2}, 'ClayMine', cm);
    
    % Strati in Faulted Section
    mySect = FaultedSection(footwall, hangingwall, faultDip, ...
                            'maxPerm', maxPerm);
    
    % Get material property distributions
    mySect = mySect.getMatPropDistr();
    
    % Get along-strike segmentation
    nSeg = getNSeg(mySect.Vcl, mySect.IsClayVcl, mySect.DepthFaulting);
    
    % Loop over realizations
    faultSections = cell(Nsim, 1);
    nSeg_fcn = nSeg.fcn;
    %smears = cell(Nsim, 1);
    parfor n=1:Nsim  
        % Generate fault object with properties for each realization
        myFaultSection = Fault2D(mySect, faultDip);
        myFault = Fault3D(myFaultSection, mySect);
        myFault = myFault.getSegmentationLength(U, nSeg_fcn);
        G = [];
        for k=1:numel(myFault.SegLen)
            % Get dependent variables
            myFaultSection = myFaultSection.getMaterialProperties(mySect, 'corrCoef', rho);
            myFaultSection.MatProps.thick = myFault.Thick;
            if isempty(G)
                G = makeFaultGrid(myFault.Thick, myFault.Disp, ...
                    myFault.Length, myFault.SegLen, U);
            end
            
            % Generate smear object with T, Tap, L, Lmax
            smear = Smear(mySect, myFaultSection, G, 1);
            
            % Place fault materials and assign cell-based properties
            myFaultSection = myFaultSection.placeMaterials(mySect, smear, G);
            
            % Extrude 2D section to fill current segment
            myFault = myFault.assignExtrudedVals(G, myFaultSection, k);
            
            faultSections{n}{k} = myFaultSection;
        end
        
        % Compute upscaled permeability distribution
        myFault = myFault.upscaleProps(G, U);
        
        % Save result
        faults{n, j} = myFault;
        %smears{n, j} = smear;        
        if mod(n, 250) == 0
            disp([num2str(n) ' realizations out of ' num2str(Nsim), ...
                ' completed.'])
        end
    end
    faultSections_all{j} = faultSections;
    sect_all{j} = mySect;

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
plotId = selectSimId('randm', faults(:,stratId), Nsim);                % simulation index

% Visualize Strati, fault thickness of 1st realization
sect_all{stratId}.plotStrati(faults{plotId,stratId}.Thick, faultDip, unit_plot);  

% Histograms for each MatProp (all sims, we select one stratigraphic layer)
% This should plot for all realizations that contain the given id.
%layerId = 2;                                            
%plotMatPropsHist(faults_all{stratId}, smears_all{stratId}, ...
%                 sect_all{stratId}, layerId, dim) 

% MatProps correlations
%[R, P] = plotMatPropsCorr(faults_all{stratId}, sect_all{stratId}, 2, dim);

% General fault materials and perm view
%faults{plotId,stratId}.plotMaterials(faultSections_all{stratId}{plotId}{1}, ...
%                                     sect_all{stratId}, unit_plot, U) 

% Plot upscaled Poro and Perm (all sims, 3 directions)
%plotUpscaledPerm(faults(:,stratId), dim)
plotUpscaledPermData(faults(:,stratId), thickness{stratId}, ...
                     vcl{stratId}, zf(stratId, :), zmax{stratId}, ...
                     idSeq(stratId, :), dim);

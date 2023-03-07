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
idSeq = [0.2, 0.45, 0.5, 0.7] ;           % strati identifiers, smear source fraction in seq.
vclSeq = [0.205 0.27 0.355 0.5275];       % avg vcl sequence
fname = 'multimodalDistros_2457_connected';         % [] (empty) to not save data.
dim = 3;

% Mandatory Input parameters
thickness = {{[20 10 20 10 40], [5 20 5 20 5 20 5 20]}, ...
             {[10 20 10 20 30 10], [5 10 15 30 20 20]}, ...
             {[10 40 30 10 10], [20 20 20 30 10]}, ...
             {[40 35 5 20], [30 20 20 30]}};
vcl       = {{[0.2 0.4 0.1 0.5 0.1], [0.5 0.2 0.4 0.1 0.5 0.3 0.4 0.05]}, ...
             {[0.4 0.1 0.5 0 0.6 0.1], [0.6 0.05 0.5 0.1 0.4 0.1]}, ...
             {[0.4 0.25 0.5 0.1 0.6], [0.55 0.1 0.5 0.2 0.6]}, ...
             {[0.1 0.7 0.4 0.8], [0.8, 0.5, 0.2, 0.7]}};
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
U.method          = 'tpfa';     % 'tpfa' in 3D
U.coarseDims      = [1 1 1];
U.flexible      = true;

% Prepare loop
nStrat = numel(thickness);
Nsim = 2000;
faults = cell(Nsim, nStrat);
faultSections_all = cell(1, nStrat);
%smears_all = cell(nStrat, 1);
sect_all = cell(nStrat, 1);
C = zeros(Nsim,nStrat,3);
Cx = nan(Nsim, 1); Cz = nan(Nsim, 1); Cy = nan(Nsim, 1);
tic
for j=1:nStrat          % For each stratigraphy
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
    
    % Get along-strike segmentation
    nSeg = getNSeg(mySect.Vcl, mySect.IsClayVcl, mySect.DepthFaulting);
    
    % Loop over realizations
    %smears = cell(Nsim, 1);
    faultSections = cell(Nsim, 1);
    nSeg_fcn = nSeg.fcn;
    parfor n=1:Nsim  
        % Generate fault object with properties for each realization
        myFaultSection = Fault2D(mySect, faultDip);
        myFault = Fault3D(myFaultSection, mySect);
        [myFault, Us{n}] = myFault.getSegmentationLength(U, nSeg_fcn);
        G = [];
        for k=1:numel(myFault.SegLen)
            % Get dependent variables
            myFaultSection = myFaultSection.getMaterialProperties(mySect, 'corrCoef', rho);
            myFaultSection.MatProps.thick = myFault.Thick;
            if isempty(G)
                G = makeFaultGrid(myFault.Thick, myFault.Disp, ...
                    myFault.Length, myFault.SegLen, Us{n});
            end

            % Generate smear object with T, Tap, L, Lmax
            smear = Smear(mySect, myFaultSection, G, 1);
            
            % Place fault materials and assign cell-based properties
            myFaultSection = myFaultSection.placeMaterials(mySect, smear, G);
            
            % Extrude 2D section to fill current segment
            myFault = myFault.assignExtrudedVals(G, myFaultSection, k);
            
            faultSections{n}{k} = myFaultSection;
            %smears{n}{k} = smear;
        end
        
        % Compute upscaled permeability distribution
        myFault = myFault.upscaleProps(G, Us{n});
            
        % Find boundary to boundary connections
        c = findSandConn(myFault.Grid.isSmear, U.method, dim, G);
        
        % Save result
        faults{n, j} = myFault;
        Cx(n) = c.bc(1); Cz(n) = c.bc(2); Cy(n) = c.bc(3);
        if mod(n, 200) == 0
            disp([num2str(n) ' realizations out of ' num2str(Nsim), ...
                ' completed.'])
        end
    end
    faultSections_all{j} = faultSections;
    sect_all{j} = mySect;
    C(:,j,1) = Cx;  C(:,j,2) = Cz;  C(:,j,3) = Cy;
    
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
sect_all{stratId}.plotStrati(faults{plotId,stratId}.Thick, faultDip, 'm'); 

% Histograms for each MatProp (all sims, we select one stratigraphic layer)
% This should plot for all realizations that contain the given id.
%layerId = 2;                                            
%plotMatPropsHist(faults_all{stratId}, smears_all{stratId}, ...
%                 sect_all{stratId}, layerId, dim) 

% MatProps correlations
%[R, P] = plotMatPropsCorr(faults_all{stratId}, sect_all{stratId}, 2, dim);

% General fault materials and perm view
faults{plotId,stratId}.plotMaterials(faultSections_all{stratId}{plotId}{1}, ...
                                     sect_all{stratId}, 'm', U) 

% Plot upscaled Poro and Perm (all sims, 3 directions)
plotUpscaledPerm(faults(:,stratId), dim)%, 'histOnly')

% Plot connected sand pathways
plotConnected(C, vclSeq, dim)

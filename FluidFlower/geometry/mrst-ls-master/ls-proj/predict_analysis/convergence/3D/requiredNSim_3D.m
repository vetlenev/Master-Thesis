%% Required number of simulations 
%  Goal: to obtain perm distros that are representative of the full 
%        parameter uncertainty. This code incrementally adds the number of
%        realizations (N) to compare the output permeability distributions 
%        as N increases.

%% Example 0: Single stratigraphic case + analysis
%
% Use parfor instead of for when running several simulations.
% 

clear
close all

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
name = {'A', 'B', 'C', 'D', 'E'};                % strati base names.
fname = 'requiredNSim_3D_5manualStrati_20k';   % [] (empty) to not save data.

% Mandatory Input parameters
thickness = {{repelem(10, 1, 10), repelem(10, 1, 10)}, ...
             {[25 25 25 25], [25 25 25 25]}, ...
             {[50 50], [50 50]}, ...
             {[5 10 15 10 20 10 10 5 15], [20 10 20 10 30 10]}, ...
             {[20 30 30 20]; [40 20 10 10 10 10]}};
vcl       = {{repmat([0.2 0.6], 1, 5), repmat([0.5 0.3], 1, 5)}, ...
             {[0.8 0.3 0.5 0], [0.3, 0.7, 0.15, 0.6]}, ...
             {[0.5 0.1], [0.5 0.1]}, ...
             {[0.3 0.6 0.1 0.7 0.2 0.8 0.3 0.9 0.1], [0.2, 0.7, 0.25, 0.8, 0.3, 0.9]}, ...
             {[0 0.4 0.3 0.6], [0.1 0.4 0.2 0.6 0.1 0.5]}};
dip       = {[0, 0], [0, 0], [10, 20], [5, -5], [0, 0]};
faultDip  = [50, 60, 75, 85, 55];

% Optional Input parameters
nl   = [numel(vcl{1}{1}), numel(vcl{1}{2}); ...
        numel(vcl{2}{1}), numel(vcl{2}{2}); ...
        numel(vcl{3}{1}), numel(vcl{3}{2}); ...
        numel(vcl{4}{1}), numel(vcl{4}{2}); ...
        numel(vcl{5}{1}), numel(vcl{5}{2})];  % just for convenience here
zf   = {[100, 100], [500, 500], [1000, 1000], [100, 100], [2000, 2000]};    % m
zmax = {{repelem(800, 1, nl(1,1)), repelem(800, 1, nl(1,2))}, ...
        {repelem(1000, 1, nl(2,1)), repelem(1000, 1, nl(2,2))}, ...
        {repelem(1000, 1, nl(3,1)); repelem(1000, 1, nl(3,2))}, ...
        {repelem(2000, 1, nl(4,1)); repelem(2000, 1, nl(4,2))}, ...
        {repelem(3000, 1, nl(5,1)); repelem(3000, 1, nl(5,2))}};
cm = {'kao', 'sme', 'ill', 'mic', 'kao'};       % predominant clay mineral
maxPerm = 5000;                                 % cap max perm? [mD]
rho = 0.6;                                      % Corr. coeff. for multiv. distr.    

% Flow upscaling options
U.useAcceleration = 1;          % requires MEX and AMGCL setup
U.method          = 'tpfa';     % 'tpfa' recommended if useAcc. = 0
U.flexible        = true;
U.coarseDims      = [1 1 1];

% Prepare loop
nStrat = numel(vcl);
nSim = [10 100 200 500 1000 2000 5000 10^4 2*10^4];
%nSim = [10 100 200 300 400 500 750 1000 2000];
log_k_md = cell(nStrat, numel(nSim));
tic
for j=1:nStrat          % Manually for each stratigraphy
    disp(['Stratigraphy ' num2str(j) ' / ' num2str(nStrat) '.'])
    fDipIt = faultDip(j);
    
    % FW and HW
    footwall = Stratigraphy(thickness{j}{1}, vcl{j}{1}, 'Dip', dip{j}(1), ...
                            'DepthFaulting', zf{j}(1), ...
                            'DepthBurial', zmax{j}{1}, 'ClayMine', cm{j});
    hangingwall = Stratigraphy(thickness{j}{2}, vcl{j}{2}, 'Dip', dip{j}(2), ...
                               'IsHW', 1, 'NumLayersFW', footwall.NumLayers, ...
                               'DepthFaulting', zf{j}(2), ...
                               'DepthBurial', zmax{j}{2}, 'ClayMine', cm{j});
    
    % Strati in Faulted Section
    mySect = FaultedSection(footwall, hangingwall, fDipIt, ...
                            'maxPerm', maxPerm);
    
    % Get material property distributions
    mySect = mySect.getMatPropDistr();
    
    % Get along-strike segmentation
    nSeg = getNSeg(mySect.Vcl, mySect.IsClayVcl, mySect.DepthFaulting);
    
    % Loop over multiple simulation numbers
    nSeg_fcn = nSeg.fcn;
    for k=1:numel(nSim)
        nSimIt = nSim(k);
        perm = nan(prod(U.coarseDims), 3, nSimIt);
        disp(['Started simulation set ' num2str(k) ' / ' num2str(numel(nSim)) ...
              ' (' num2str(nSimIt) ' realizations)...'])
        parfor n=1:nSimIt    % loop for each realization in a given set 
            % Generate fault object with properties for each realization
            myFaultSection = Fault2D(mySect, fDipIt);
            myFault = Fault3D(myFaultSection, mySect);
            myFault = myFault.getSegmentationLength(U, nSeg_fcn);
            G = [];
            for q = 1:numel(myFault.SegLen)
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
                myFault = myFault.assignExtrudedVals(G, myFaultSection, q);
            
            end
            % Compute upscaled permeability distribution
            myFault = myFault.upscaleProps(G, U);
            
            % Save result
            perm(:, :, n) = myFault.Perm;
            
            if mod(n, 500) == 0
                disp([num2str(n) ' realizations out of ' num2str(nSimIt), ...
                      ' completed.'])
            end
        end
        log_k_md{j, k} = log10(perm ./ (milli*darcy));
    
        disp(['Simulation ' num2str(k) ' done. ' num2str(nStrat - j) ...
              ' stratigraphies remaining.'])
        disp('-----------------------------------------------------------')
    end
    
    disp(['Stratigraphy ' num2str(j) ' finished.'])
    disp('***************************************************************')
end
telapsed = toc;

% Save data?
if ~isempty(fname)
    disp(['ATTENTION: data saved in: ' pwd ' with filename ' fname])    
    save([fname '.mat']) %,'-v7.3') % larger than 2GB
end


%% Output analysis
nedg = 51;
edges = linspace(-7, 4, nedg); 
nbin = nedg - 1;

% load separate data if each strat run and saved separately
% log_k_md_e = log_k_md(end,:);
% log_k_md_a = log_k_md(1,:);
% log_k_md_b = log_k_md(2,:);
% log_k_md_c = log_k_md(3,:);
% log_k_md_d = log_k_md(4,:);
% log_k_md(1,:) = log_k_md_a;
% log_k_md(2,:) = log_k_md_b;
% log_k_md(3,:) = log_k_md_c;
% log_k_md(4,:) = log_k_md_d;
% log_k_md(5,:) = log_k_md_e;

% Compare perm distros
maex = nan(nStrat, numel(nSim)-1); maey = maex; maez = maex;
ks2x = nan(nStrat, numel(nSim)-1); ks2y = ks2x; ks2z = ks2x;
for j = 1:nStrat
    log_kxx_last = reshape(log_k_md{j, end}(:,1,:), prod(U.coarseDims)*nSim(end), 1);
    log_kyy_last = reshape(log_k_md{j, end}(:,2,:), prod(U.coarseDims)*nSim(end), 1);
    log_kzz_last = reshape(log_k_md{j, end}(:,3,:), prod(U.coarseDims)*nSim(end), 1);
    p = [histcounts(log_kxx_last, edges, 'normalization', 'probability')', ...
         histcounts(log_kyy_last, edges, 'normalization', 'probability')', ...
         histcounts(log_kzz_last, edges, 'normalization', 'probability')'];
    idx = find(p(:,1)); idy = find(p(:,2)); idz = find(p(:,3)); % only count bins where p > 0
    for k = 1:(numel(nSim) - 1)
        log_kxx = reshape(log_k_md{j, k}(:,1,:), prod(U.coarseDims)*nSim(k), 1);
        log_kyy = reshape(log_k_md{j, k}(:,2,:), prod(U.coarseDims)*nSim(k), 1);
        log_kzz = reshape(log_k_md{j, k}(:,3,:), prod(U.coarseDims)*nSim(k), 1);
        % MAE
        pit = [histcounts(log_kxx, edges, 'normalization', 'probability')', ...
               histcounts(log_kyy, edges, 'normalization', 'probability')', ...
               histcounts(log_kzz, edges, 'normalization', 'probability')'];
        maex(j, k) = sum(abs(pit(idx, 1) - p(idx, 1))) / numel(idx);
        maey(j, k) = sum(abs(pit(idy, 2) - p(idy, 2))) / numel(idy);
        maez(j, k) = sum(abs(pit(idz, 3) - p(idz, 3))) / numel(idz);
        
        % KS2
        ks2x(j, k) = kstest2(log_kxx, log_kxx_last);
        ks2y(j, k) = kstest2(log_kyy, log_kyy_last);
        ks2z(j, k) = kstest2(log_kzz, log_kzz_last);
    end
end

% Plot utils
sz = [14, 12];
latx = {'Interpreter','latex'};

% Plot
colrs = [0.35 0.35 0.35; 1 0 0; 0 0 1];
plotId = [1 3 5 8 9];
%plotId = 1:3; 
fh = figure(1);
tiledlayout(numel(plotId), nStrat, 'Padding', 'compact', 'TileSpacing', 'compact');
tileIds = reshape(1:nStrat*numel(plotId), nStrat, numel(plotId))';
for j=1:nStrat
    log_kxx_last = reshape(log_k_md{j, end}(:,1,:), prod(U.coarseDims)*nSim(end), 1);
    log_kyy_last = reshape(log_k_md{j, end}(:,2,:), prod(U.coarseDims)*nSim(end), 1);
    log_kzz_last = reshape(log_k_md{j, end}(:,3,:), prod(U.coarseDims)*nSim(end), 1);
    %limx = [floor(log10(min(min(k_md{j, end})))), ...
    %        ceil(log10(max(max(k_md{j, end}))))];
    %edges = logspace(limx(1), limx(2), nbins);
    probs = [histcounts(log_kxx_last, edges); ...
             histcounts(log_kyy_last, edges); ...
             histcounts(log_kzz_last, edges)]' ./ nSim(end);
    limy = floor(log10(min(probs(probs > 0))));
    for k=1:numel(plotId)
        log_kxx = reshape(log_k_md{j, plotId(k)}(:,1,:), [], 1);
        log_kyy = reshape(log_k_md{j, plotId(k)}(:,2,:), [], 1);
        log_kzz = reshape(log_k_md{j, plotId(k)}(:,3,:), [], 1);
        nexttile(tileIds(k, j))
        hold on
        histogram(log_kxx, edges, 'Normalization', 'probability', ...
            'FaceColor', colrs(1, :), 'EdgeColor', colrs(1,:), ...
            'FaceAlpha', 1, 'DisplayName', '$k_{xx}$')
        histogram(log_kyy, edges, 'Normalization', 'probability', ...
            'FaceColor', colrs(2, :), 'EdgeColor', 'none', ...
            'FaceAlpha', .7, 'DisplayName', '$k_{yy}$')
        histogram(log_kzz, edges, 'Normalization', 'probability', ...
            'FaceColor', colrs(3, :), 'EdgeColor', 'none', ...
            'FaceAlpha', .5, 'DisplayName', '$k_{zz}$')
        if k == 1
            xlabel('$\log_{10}(k$ [mD])', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            if j == 1
                title([name{j} ', $N_\mathrm{sim}$ = ' num2str(nSim(plotId(k)))], latx{:}, ...
                   'fontSize', sz(2))
                leg = legend(latx{:}, 'fontSize', sz(2), 'location', 'northwest');
                set(leg.BoxFace, 'ColorType','truecoloralpha', ...
                    'ColorData', uint8(255*[1;1;1;.6]));
            else
                title(name{j}, latx{:}, 'fontSize', sz(2))
            end
        elseif j == 1
            if nSim(plotId(k)) < 10000
                text(-3.5, 0.5, num2str(nSim(plotId(k))), latx{:}, 'fontSize', 11)
            else
                text(-3.5, 0.5, num2str(nSim(plotId(k)), 3), latx{:}, 'fontSize', 11)
            end
        end
        %xlim([10^limx(1, 1) 10^limx(1, 2)])
        %xticks(10.^(limx(1, 1):2:limx(1, 2)));
        xlim([edges(1) edges(end)])
        xticks(-6:2:4)
        %ylim([10^limy 1])
        %yticks(10.^(limy:0))
        if j == 1, ylim([0 0.65]); yticks(0:.1:.6)
        elseif j == 2, ylim([0 0.4]); yticks(0:.1:.4)
        elseif j == 3, ylim([0 0.85]); yticks(0:.2:1)   
        elseif j == 4, ylim([0 0.6]); yticks(0:.1:.6)
        else, ylim([0 0.5]); yticks(0:.1:.6)
        end
        grid on
        %set(gca,'XScale','log','YScale','log')
        %set(gca, 'XScale', 'log');
    end    
end
set(fh, 'position', [200, 200, 850, 600]);
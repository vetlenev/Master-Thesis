function [G, rock, myFault, mySect, faultSections, Us, CG] = ...
                                            upscaleFaultPerm(opt, makeplots)
%
%
%
if strcmp(opt.fault,  'test')
    G = computeGeometry(cartGrid([100, 100], [5, 100]));
    
    kvals = [1e-3 100]*(milli*darcy);
    rock.perm = repelem(kvals(2), G.cells.num, 1);
    id_lowperm = rand(G.cells.num, 1) <= 0.5;
    rock.perm(id_lowperm) = kvals(1);
    
    rock.poro = repelem(0.3, G.cells.num, 1);
    rock.poro(id_lowperm) = 0.15;
    myFault.Grid.perm = rock.perm;
    p  = partitionCartGrid(G.cartDims, [1 1]);
    CG = generateCoarseGrid(G, p);
    myFault.Perm = [0 0];
    myFault.Perm(end) = myUpscalePermDim(G, CG, rock, 'dim', 2);
    mySect = [];
    faultSections = [];
    Us = [];
    
elseif strcmp(opt.fault, 'predict')
    rho               = 0.6;        % Corr. coeff. for multivariate distributions
    U.useAcceleration = 1;          % 1 requires MEX setup, 0 otherwise (slower for MPFA).
    U.method          = 'mpfa';     % 'tpfa' recommended if useAcceleration = 0
    U.outflux         = 0;          % compare outflux of fine and upscaled model
    U.ARcheck         = 0;          % check if Perm obtained with grid with aspect ratio of
                                    % only 5 gives same output.
    %Nsim              = 1000;       % Number of simulations/realizations
    
    footwall = Stratigraphy(opt.thick{1}, opt.vcl{1}, 'Dip', opt.dip(1), ...
                           'DepthFaulting', opt.zf(1), 'DepthBurial', opt.zmax{1});
    hangingwall = Stratigraphy(opt.thick{2}, opt.vcl{2}, 'Dip', opt.dip(2), ...
                          'IsHW', 1, 'NumLayersFW', footwall.NumLayers, ...
                          'DepthFaulting', opt.zf(2), 'DepthBurial', opt.zmax{2});
    
    % Instantiate FaultedSection object (Strati in Faulted Section)
    mySect = FaultedSection(footwall, hangingwall, opt.fDip, 'maxPerm', opt.maxPerm);
    
    % Get material distributions
    mySect = mySect.getMatPropDistr();
    
    % Generate intermediate variable samples, calculate smear dimensions 
    % and upscale permeability
    %faults = cell(Nsim, 1);
    %smears = cell(Nsim, 1);
    %tstart = tic;
    %parfor n=1:Nsim    % parfor allowed if you have the parallel computing toolbox
        myFault = Fault2D(mySect, opt.fDip);
        
        % Get material property (intermediate variable) samples
        myFault = myFault.getMaterialProperties(mySect, 'corrCoef', rho);
        
        % Generate smear object with T, Tap, L, Lmax
        smear = Smear(mySect, myFault, 1);
        
        % Compute upscaled permeability distribution
        [myFault, G] = myFault.upscaleSmearPerm(mySect, smear, U);
        
        % Save result
        %faults{n} = myFault;
        %smears{n} = smear;
        %if mod(n, 100) == 0
        %    disp(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
        %end
    %end
    %telapsed = toc(tstart);
    rock.poro = myFault.Grid.poro;
    rock.perm = myFault.Grid.perm;
    
elseif strcmp(opt.fault, 'predict_3D')
    rho               = 0.6;        % Corr. coeff. for multivariate distributions
    U.useAcceleration = 1;          % 1 requires MEX setup, 0 otherwise (slower for MPFA).
    U.method          = 'tpfa';     % 'tpfa' recommended if useAcceleration = 0
    U.coarseDims      = [1 1 1];          
    U.flexible        = true; 
    
    footwall = Stratigraphy(opt.thick{1}, opt.vcl{1}, 'Dip', opt.dip(1), ...
                           'DepthFaulting', opt.zf(1), 'DepthBurial', opt.zmax{1});
    hangingwall = Stratigraphy(opt.thick{2}, opt.vcl{2}, 'Dip', opt.dip(2), ...
                          'IsHW', 1, 'NumLayersFW', footwall.NumLayers, ...
                          'DepthFaulting', opt.zf(2), 'DepthBurial', opt.zmax{2});
    
    % Instantiate FaultedSection object (Strati in Faulted Section)
    mySect = FaultedSection(footwall, hangingwall, opt.fDip, 'maxPerm', opt.maxPerm);
    % Get material property distributions
    mySect = mySect.getMatPropDistr();
    % Get strike-parallel segmentation
    nSeg = getNSeg(mySect.Vcl, mySect.IsClayVcl, mySect.DepthFaulting);
    
    % Instantiate fault section and get segmentation for this realization
    myFaultSection = Fault2D(mySect, opt.fDip);
    myFault = Fault3D(myFaultSection, mySect);
    if U.flexible
        [myFault, Us] = myFault.getSegmentationLength(U, nSeg.fcn);
    else
        myFault = myFault.getSegmentationLength(U, nSeg.fcn);
    end
    G = [];
    faultSections = cell(numel(myFault.SegLen), 1);
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
        faultSections{k} = myFaultSection;
        %smears{k} = smear;
    end
    
    % Compute 3D upscaled permeability distribution
    if U.flexible
        [myFault, CG] = myFault.upscaleProps(G, Us);
    else
        [myFault, CG] = myFault.upscaleProps(G, U);
    end
    % Assign rock props
    rock.poro = myFault.Grid.poro;
    rock.perm = myFault.Grid.perm;
end


%% Plots

% Grid and perm
if makeplots
    if strcmp(opt.fault, 'predict') 
        mySect.plotStrati(myFault.MatProps.thick, opt.fDip, 'm'); 
    elseif strcmp(opt.fault, 'predict_3D')
        mySect.plotStrati(myFault.Thick, opt.fDip, 'm'); 
    end
    
    f=figure(2); colormap(copper);
    plotCellData(G, log10(rock.perm(:,1)/(milli*darcy)),'edgecolor', 'none'), colorbar
    %axis equal
    f.Position = [100 100 100 500];
    if strcmp(opt.fault, 'predict_3D')
        ax = gca;
        ax.DataAspectRatio = [0.05 1 1];
        ax.ZDir = 'normal';
        view([30 20])
        f.Position = [100 100 500 500];
        
        if U.flexible
            myFault.plotMaterials(faultSections{1}, mySect, 'm', Us)
        else
            myFault.plotMaterials(faultSections{1}, mySect, 'm', U)
        end
    end
end

end
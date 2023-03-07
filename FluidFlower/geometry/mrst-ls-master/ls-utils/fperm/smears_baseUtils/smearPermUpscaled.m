function [K, flux, V, D] = smearPermUpscaled(G, permG, f, plotFigs, ...
                                             gcellsCheck, kbc, kmethod, kacc, outfl, outUps)
%
% Perm is the permeability tensor (2x2, x and z dimensions) as obtained
% from a single-phase, incompressible flow simulation to steady state
% using a Mimetic solver (consistent also for grids that are not
% k-orthogonal) and periodic boundary conditions, to get off-diagonal
% terms.
%

%mrstModule add mpfa mimetic coarsegrid upscaling incomp; 

L = max(G.faces.centroids) - min(G.faces.centroids);
fluid = initSingleFluid('mu' , 1, 'rho', 1);
%rock.poro = poroG;
rock.perm = permG;

switch kbc          % boundary conditions to determine upscaled permeability
    case 'sealed'   % p drop in each axial dim with other boundaries sealed
        
        % compute equivalent/upscaled perm according to input method
        if strcmp(kmethod, 'tpfa')
            p2 = partitionCartGrid(G.cartDims, [1 1]);
            CG2 = generateCoarseGrid(G, p2);
            K = diag(upscalePerm(G, CG2, rock, 'method', kmethod));
            
        elseif strcmp(kmethod, 'mpfa')
            Dp{1} = 5*barsa;
            if kacc == 1
                inB = 'mex';
            else
                inB = [];
            end
            hTmp = computeMultiPointTrans(G, rock, 'invertBlocks', inB);
            psolver = @(state0, G, fluid, bc) incompMPFA(state0, G, hTmp, fluid, 'bc', bc);
            K = diag(myupscalePermeabilityFixed(G, Dp{1}, psolver, fluid, L));
            
        elseif strcmp(kmethod, 'mimetic')
            Dp{1} = 5*barsa;
            Sp2   = computeMimeticIP(G, rock);
            psolver = @(state0, G, fluid, bc) incompMimetic(state0, G, Sp2, ...
                                                            fluid, 'bc', bc);
            K = diag(myupscalePermeabilityFixed(G, Dp{1}, psolver, fluid, L));
        end
        
        % compare outflux if requested
        if outfl == 1
            % Pressure drop along the axial (Y) dimension, with sealed lateral
            % boundaries. Remember that MRST uses SI units (perm m^2)
            
            % Half-transmissibilities
            hT = computeTrans(G, rock);
            
            % Coarse grid
            upscaled = outUps;
            disp('________________________________________________________________')
            disp(' Pressure drop with sealed boundaries: ')
            disp(['Coarse grid cells along axial direction is ' num2str(upscaled(2))]);
            G_ups = computeGeometry(cartGrid(upscaled, [f.T f.D]));
            p = partitionCartGrid(G.cartDims, upscaled);
            CG = generateCoarseGrid(G, p);
            
            % Define coarse grid permeability
            crock.perm = repmat([K(1) K(4)], CG.cells.num, 1);
            
            % Define fluid
            fluid = initSingleFluid('mu', 1*centi*poise, ...
                'rho', 1014*kilogram/meter^3);                      % water
            
            % Define initial conditions
            watCol = 50*meter;
            topy   = 1000*meter;
            g      = 9.861;
            p_r    = (watCol + topy + min(G.cells.centroids(:,2)))*fluid.rhoWS*g;
            [z_0, z_max] = deal(min(G.cells.centroids(:,2)), ...
                max(G.cells.centroids(:,2)));
            equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0, z_max], p_r);
            p0     = flipud(reshape(deval(equil, G.cells.centroids(:,2)), [], 1));
            state0 = initResSol(G, p0, 1); clear equil
            
            p_s = (watCol + topy + min(G_ups.cells.centroids(:,2)))*fluid.rhoWS*g;
            [z_0u, z_maxu] = deal(min(G_ups.cells.centroids(:,2)), ...
                max(G_ups.cells.centroids(:,2)));
            if z_0u ~= z_maxu
                equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0u, z_maxu], p_s);
                p0u = flipud(reshape(deval(equil, G_ups.cells.centroids(:,2)), [], 1));
                clear equil
            elseif G_ups.cells.num == 1
                p0u = (watCol + topy + G_ups.cells.centroids(:,2))*fluid.rhoWS*g;
            end
            state0c = initResSol(G_ups, p0u, 1);
            
            % Fine-scale problem
            bc        = pside([], G, 'North', 0.75*min(p0));
            faces     = bc.face;
            bc        = pside(bc, G, 'South',  1.25*max(p0));
            x         = incompTPFA(state0, G, hT, fluid, 'bc', bc);
            if strcmp(kmethod, 'mpfa')
                xm   = incompMPFA(state0, G, hTmp, fluid, 'bc', bc);
                str  = 'Sum outflux on fine scale (MPFA)  : ';
            else
                S    = computeMimeticIP(G, rock);
                xm   = incompMimetic(state0, G, S, fluid, 'bc', bc);
                str  = 'Sum outflux on fine scale (Mimetic)  : ';
            end
            
            % Coarse-scale problem
            bc_ups    = pside([], G_ups, 'North', 0.75*min(p0));
            faces_ups = bc_ups.face;
            bc_ups    = pside(bc_ups, G_ups, 'South', 1.25*max(p0));       
            T_ups     = computeTrans(G_ups, crock);
            x_ups     = incompTPFA(state0c, G_ups, T_ups, fluid, 'bc', bc_ups);
            
            % Results
            flux1 = sum(x.flux(faces));
            flux2 = sum(xm.flux(faces));
            flux3 = sum(x_ups.flux(faces_ups));
            disp(['Sum outflux on fine scale (TPFA)  : ', num2str(flux1)]);
            disp([str, num2str(flux2)]);
            disp(['Sum outflux on coarse scale (TPFA with p drop, ' ...
                  num2str(kmethod) ' k_ups) : ', num2str(flux3)]);
            disp('___________________________________________________________________');
            flux.axialDropType = {'Fine scale TPFA', 'Fine scale Mimetic', ...
                                'Coarse scale TPFA with p drop, mimetic k'};
            flux.axialDropVals = [flux1, flux2, flux3];
            
        end
        
    case 'open'     % periodic boundary conditions
        % Structures with boundary conditions
        d = G.griddim;
        [bcl,bcr, Dp]=deal(cell(d,1));
        bcsides = {'XMin', 'XMax'; 'YMin', 'YMax'; 'ZMin', 'ZMax'};
        for j = 1:d
            bcl{j} = pside([], G, bcsides{j, 1}, 0);
            bcr{j} = pside([], G, bcsides{j, 2}, 0);
            Dp{j}  = 0;
        end
        Dp{1} = 5*barsa;
        
        % Periodic grid
        [Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, Dp);
        
        % Compute upscaled permeability
        Sp2     = computeMimeticIPGp(G, Gp, rock);
        psolver = @(state0, Gp, fluid, bcp) incompMimetic(state0, Gp, Sp2, ...
                                                          fluid, 'bcp', bcp);
        K = myupscalePermeabilityPeriodic(Gp, bcp, Dp{1}, psolver, fluid, L);    % m^2
        [V,D] = eig(K);
        assert(all(eig(K)>0))
        
        if outfl == 1 
            % P drop, open boundaries at constant p
            % Half-transmissibilities
            hT = computeTrans(G, rock);
            
            % Coarse grid
            upscaled = outUps;
            disp('________________________________________________________________')
            disp(' Pressure drop in axial (y) direction, open lateral boundaries: ')
            disp(['Coarse grid cells along axial direction is ' num2str(upscaled(2))]);
            G_ups = computeGeometry(cartGrid(upscaled, [f.T f.D]));
            p = partitionCartGrid(G.cartDims, upscaled);
            CG = generateCoarseGrid(G, p);
            
            % Define coarse grid permeability
            crock.perm = repmat([K(1) K(2) K(4)], CG.cells.num, 1);
            
            % Define fluid
            fluid = initSingleFluid('mu', 1*centi*poise, ...
                'rho', 1014*kilogram/meter^3);                      % water
            
            % Define initial conditions
            watCol = 50*meter;
            topy   = 1000*meter;
            g      = 9.861;
            p_r    = (watCol + topy + min(G.cells.centroids(:,2)))*fluid.rhoWS*g;
            [z_0, z_max] = deal(min(G.cells.centroids(:,2)), ...
                max(G.cells.centroids(:,2)));
            equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0, z_max], p_r);
            p0     = flipud(reshape(deval(equil, G.cells.centroids(:,2)), [], 1));
            state0 = initResSol(G, p0, 1); clear equil
            
            p_s = (watCol + topy + min(G_ups.cells.centroids(:,2)))*fluid.rhoWS*g;
            [z_0u, z_maxu] = deal(min(G_ups.cells.centroids(:,2)), ...
                max(G_ups.cells.centroids(:,2)));
            if z_0u ~= z_maxu
                equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0u, z_maxu], p_s);
                p0u = flipud(reshape(deval(equil, G_ups.cells.centroids(:,2)), [], 1));
                clear equil
            elseif G_ups.cells.num == 1
                p0u = (watCol + topy + G_ups.cells.centroids(:,2))*fluid.rhoWS*g;
            end
            state0c = initResSol(G_ups, p0u, 1);
            
            % Fine scale problem
            pTop = 0;
            pBot = 10*max(p0);
            bfac = boundaryFaces(G);
            fSide = bfac(G.faces.centroids(bfac, 1) < .001);
            fSideCentr = G.faces.centroids(fSide, :);
            
            p_s    = (watCol + topy + min(fSideCentr(:,2)))*fluid.rhoWS*g;
            [z_0, z_max] = deal(min(fSideCentr(:,2)), max(fSideCentr(:,2)));
            equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0, z_max], p_s);
            pSideVals = flipud(reshape(deval(equil, fSideCentr(:,2)), [], 1));  clear equil
            
            
            bc = pside([], G, 'North', pTop);
            faces = bc.face;
            bc = pside(bc, G, 'South', pBot);
            bc = pside(bc, G, 'East', pSideVals);
            bc = pside(bc, G, 'West', pSideVals);
            
            x  = incompTPFA(state0, G, hT, fluid, 'bc', bc);
            S  = computeMimeticIP(G, rock);
            xm = incompMimetic(state0, G, S, fluid, 'bc', bc);
            
            % Coarse-scale problem
            bfacu = boundaryFaces(G_ups);
            fSideu = bfacu(G_ups.faces.centroids(bfacu, 1) < .001);
            fSideCentru = G_ups.faces.centroids(fSideu, :);
            
            % hydrostatic pressure on open lateral face
            p_su    = (watCol + topy + min(fSideCentru(:,2)))*fluid.rhoWS*g;
            [z_0u, z_maxu] = deal(min(fSideCentru(:,2)), max(fSideCentru(:,2)));
            if z_0u ~= z_maxu
                equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0u, z_maxu], p_su);
                pSideValsu = flipud(reshape(deval(equil, fSideCentru(:,2)), [], 1));  clear equil
            elseif G_ups.cells.num == 1
                pSideValsu = (watCol + topy + fSideCentru(:,2))*fluid.rhoWS*g;
            end
            
            % Apply boundary conditions
            bc_ups = pside([], G_ups, 'North', pTop);
            faces_ups = bc_ups.face;
            bc_ups = pside(bc_ups, G_ups, 'South', pBot);
            bc_ups = pside(bc_ups, G_ups, 'East', pSideValsu);
            bc_ups = pside(bc_ups, G_ups, 'West', pSideValsu);
            
            % Solve
            T_ups    = computeTrans(G_ups, crock);
            x_ups    = incompTPFA(state0c, G_ups, T_ups, fluid, 'bc', bc_ups);
            S_ups    = computeMimeticIP(G_ups, crock);
            x_upsm   = incompMimetic(state0c, G_ups, S_ups, fluid, 'bc', bc_ups);
            
            % Results
            flux1 = sum(x.flux(faces));
            flux2 = sum(xm.flux(faces));
            flux3 = sum(x_ups.flux(faces_ups));
            flux4 = sum(x_upsm.flux(faces_ups));
            disp(['Sum outflux on fine scale (TPFA): ', num2str(flux1)]);
            disp(['Sum outflux on fine scale (Mimetic): ', num2str(flux2)]);
            disp(['Sum outflux on coarse scale (TPFA with periodic mimetic k): ', num2str(flux3)]);
            disp(['Sum outflux on coarse scale (Mimetic with periodic mimetic k): ', num2str(flux4)]);
            
            flux.axialDropOpenLatType = {'Fine scale TPFA', 'Fine scale mimetic', ...
                                          'Coarse scale TPFA with periodic, mimetic k', ...
                                          'Coarse scale mimetic with periodic, mimetic k'};
            flux.axialDropOpenLatVals = [flux1, flux2, flux3 flux4];
        end
end

if gcellsCheck == 1 % check if same perm with smaller aspect ratio. This can be slow.
   assert(G.griddim == 2)               % way too many elements in 3D
   assert(strcmp(kbc, 'sealed'));       % comparison for sealed BCs
   ar = G.xzFaceDim(2)/G.xzFaceDim(1); 
   ar2 = 3;
   if ar > ar2
       disp('Checking if cell aspect ratio is too high for flow-based upscaling')
       Lz = G.cartDims(2)*G.xzFaceDim(2);
       nels = [G.cartDims(1) round(Lz/(ar2*G.xzFaceDim(1)))];
       G2 = computeGeometry(cartGrid([nels(1), nels(2)], [f.T f.D]));
       disp(['Finer grid with cell aspect ratio = ' num2str(ar2) ' created.'])
       disp(['Finer grid has a total of ' num2str(G2.cells.num) ' cells.' ...
             ' Now computing cell distances...'])
        
       % assign perm based on cell centroid distance to closest cell
       distx = pdist2(G2.cells.centroids(:,1), G.cells.centroids(:,1));
       distz = pdist2(G2.cells.centroids(:,end), G.cells.centroids(:,end));
       dist = sqrt(distx.^2 + distz.^2);
       disp('Distances computed. Now calculating MPFA perm with finer z mesh...') 
        
       % Compute 2nd perm
       [~, mapG2toG] = min(dist,[],2);
       rock2.perm = rock.perm(mapG2toG, :);
       hTmp2 = computeMultiPointTrans(G2, rock2, 'invertBlocks', 'mex');
       psolver = @(state0, G, fluid, bc) incompMPFA(state0, G, hTmp2, fluid, 'bc', bc);
       fluid = initSingleFluid('mu' , 1, 'rho', 1);
       K2 = diag(myupscalePermeabilityFixed(G2, Dp{1}, psolver, fluid, L));
       disp(['Finer z mesh perm [x, z], in m^2, is [' num2str(K2(1)) ' ' num2str(K2(end)) ']'])
       disp(['Ratio Perm(finerGrid) / Perm(originalGrid) [x, z] is [' ...
             num2str(K2(1)/K(1)) ' ' num2str(K2(end)/K(end)) ']'])
       
       % Plot
       figure(15); plotCellData(G2, log10(rock2.perm(:,end)./(milli*darcy)), 'facea', 1, 'edgea', 0); 
       colormap(copper); title('Finer Z Grid vertical Perm')
       xlim([0 f.T]); ylim([0 f.D]); c = colorbar;
   else
       disp(['aspect ratio < ' num2str(ar2) '. Check not needed'])
   end
    
end


if outfl == 1 && plotFigs == 1
    figure(30); subplot(1,3,1)
    plotCellData(G,log10(rock.perm(:,1)/(milli*darcy)),'EdgeColor','none');
    colorbar; axis equal; xlim([0 f.T]); ylim([0 f.D]);
    
    subplot(1,3,2)
    plotCellData(G, p0,'EdgeColor','none'); c = colorbar;
    axis equal, xlabel('x [m]'), ylabel('y [m]');
    xlim([0 f.T]); ylim([0 f.D]);
    
    if strcmp(kbc, 'open')
        subplot(1,3,3)
        plot(pSideVals/1e+6, G.faces.centroids(fSide, 2), '.-k')
        hold on
        plot(pSideValsu/1e+6, G_ups.faces.centroids(fSideu, 2), 'sb')
        xlabel('p applied to laterals [MPa]'), ylabel('y coordinate [m]');
        legend('fine grid', 'coarse grid'); grid on
        ylim([0 f.D]);
    end
    
    streamColor = [69, 118, 255]/255;
    figure(31)
    h1 = subplot(1,2,1);
    plotCellData(G, log10(rock.perm(:,end)./(milli*darcy)), 'facea', 1, 'edgea', 0);
    c = colorbar; c.Label.String = 'log_{10} k [mD]';
    set(h1, 'colormap', [0 0 0; 252, 244, 189]/255)
    hold on;
    %cells = (G.cartDims(1):G.cartDims(1)-1:prod(G.cartDims)-1)';
    startCells = (1:10:G.cartDims(1))';
    s1 = streamline(pollock(G, x, startCells, 'substeps', 1));
    x.flux = -x.flux;
    s2 = streamline(pollock(G, x, startCells, 'substeps', 1));
    set([s1; s2],'Color',streamColor);
    hold off;
    xlim([0 f.T]); ylim([0 f.D]);
    title('Streamlines over fault k')
    
    h2 = subplot(1,2,2);
    plotCellData(G, xm.pressure/1e+6, 'facea', 1, 'edgea', 0);
    c = colorbar; c.Label.String = 'p [MPa]';
    set(h2, 'colormap', flipud(hot))
    hold on;
    s1 = streamline(pollock(G, xm, startCells, 'substeps', 1));
    xm.flux = -xm.flux;
    s2 = streamline(pollock(G, xm, startCells, 'substeps', 1));
    set([s1; s2],'Color',streamColor);
    hold off;
    xlim([0 f.T]); ylim([0 f.D]);
    title('Streamlines over p field')
    
%elseif plotFigs == 1
%    flux = [];
%    figure(30); subplot(1,2,1)
%    plotCellData(G,log10(rock.perm(:,1)/(milli*darcy)),'EdgeColor','none');
%    colorbar; axis equal; xlim([0 f.T]); ylim([0 f.D]);
end


end
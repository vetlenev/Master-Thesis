function [G, rock, ucids, fault] = getRockParams_bo(mesh, G, opt, fig_perm)                                 
%
% NOTE: permeabilities and porosities in cells outside the reservoir are
%       also anisotropic and given according to their lithology. They
%       remain unchanged regardless of petrophysics properties considered
%       in the storage reservoir and fault.
%       Because some cases introduce anisotropic permeabilities in the
%       reservoir (aligned with global x,y,z, 3 values), and anisotropic 
%       permeabilities for the fault (whose principal directions are not 
%       aligned with global x,y,z, so 6 values), 6 values are given for
%       each cell. 
%       · For isotropic, homogeneous layers:
%         K11=Kxx=K22=Kyy=K33=Kzz (rock.perm col. 1,4,6), K12=K13=K23=0
%       · For anisotropic layers aligned with global x,y,z:
%         K11=K22 \neq 0, K33 \neq 0, K12=K13=K23=0
%       · For anisotropic "layers" (fault), whose input permeabilities,
%         always given as the principal K values, are NOT aligned with
%         global y,z (permeability along strike is already aligned with 
%         global x, the extruded axis)
%         K_rot = R'*[K11 0 0; 0 K22 0; 0 0 K33]*R 
%         where R = [1 0 0; 0 c(90-fault_dip) s; 0 -s c]
%         rock.perm col. 1 (K11_glob=Kxx) = K_rot(1,1) = K11,
%         rock.perm col. 2 (Kxy) = col. 3 (Kxz) = 0; col 5 (Kyz) \neq 0
%         rock.perm col. 4 (Kyy, across fault) \neq K22
%         rock.perm col. 6 (Kzz, updip) \neq K33
%

%% Parameter assignment
% base porosity and permeability. See permTensor.m
rock.poro = repmat(.1, [G.cells.num, 1]);
rock.perm = repmat(300*milli*darcy(), [G.cells.num, 6]);
rock.perm(:,[2 3 5]) = 0;

% Get cell ids of each unit
obtained.ucids    = 1;
obtained.fnodeset = 1;
tlim.bot          = 181.2952;             % throw at fault (thick) bottom
obtained.elemOnTop = 1;
if obtained.ucids == 1
    switch mesh.type
%         case 'coarse'
%             %ucids = load('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/extr_tri/ucids.mat');
%             ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/ucids_WideX.mat');
%         case 'coarseUpper'
%             ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/ucids_wideX_coarseUpper.mat');
%             ucids.upper = 1;
%         case 'ref300'
%             ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref300/ucids_ref300.mat');
%         case 'ref200'
%             ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref200/ucids_ref200.mat');
%         case 'ref100'
%             ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref100/ucids_ref100.mat');
%         case 'ref50'
%             ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref50/ucids_ref50.mat');
%         case 'ref25'
%             ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref25/ucids_ref25_upper.mat');
%             ucids.upper = 1;
%         case 'ref12.5'
%             ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref12.5/ucids_ref12.5_upper.mat');
%             ucids.upper = 2;
%         case 'sc21'
%             ucids = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/sc2/init/ucids_sc21_upper.mat');
%             ucids.upper = 1;
%             ucids.unit_cell_ids{5} = sort(ucids.unit_cell_ids{5});
%             ucids.unit_cell_ids{6} = sort(ucids.unit_cell_ids{6});
%             ucids.unit_cell_ids{7} = sort(ucids.unit_cell_ids{7});
%             ucids.unit_cell_ids{8} = sort(ucids.unit_cell_ids{8});
        case 'sc2'
            if ischar(opt.studyCase) && strcmp(opt.studyCase, 'sc2_predict_theta30_coarse') ...
               || strcmp(opt.studyCase, 'sc2_predict_theta30_coarse_nohyst') 
                ucids = load(fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/2.5Dmesh/extr_tri/sc2/throwDiv/29layers/ucids_sc2_upper.mat'));
            else
                ucids = load(fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/2.5Dmesh/extr_tri/sc2/throwDiv/ucids_sc2_upper.mat'));
            end
            ucids.upper = 1;
            ucids.unit_cell_ids{25} = sort(ucids.unit_cell_ids{25});
            ucids.unit_cell_ids{26} = sort(ucids.unit_cell_ids{26});
            ucids.unit_cell_ids{27} = sort(ucids.unit_cell_ids{27});
            ucids.unit_cell_ids{28} = sort(ucids.unit_cell_ids{28});
            ucids.unit_cell_ids{29} = sort(ucids.unit_cell_ids{29});
            ucids.unit_cell_ids{30} = sort(ucids.unit_cell_ids{30});
            ucids.unit_cell_ids{31} = sort(ucids.unit_cell_ids{31});
            ucids.unit_cell_ids{32} = sort(ucids.unit_cell_ids{32});
            ucids.unit_cell_ids{33} = sort(ucids.unit_cell_ids{33});
    end
    ucids.unit_cell_ids{end} = sort(ucids.unit_cell_ids{end});  % IMPORTANT: sort cell id for the fault for region evaluation in simulateScheduleAD
    
    % Layer poro and perm (top to bottom, except injection reservoir & fault)
    [rock, ucids, idSealLyrs, idResLyrs] = setLayerPoroPerm(mesh, rock, ucids, opt);
    
    % Fault poro and perm
    [rock, fault] = setFaultPoroPerm(mesh, rock, ucids, opt, G, tlim, idSealLyrs, idResLyrs, fig_perm);
    
    % Update top-layer centroids (triangle Grid) in injection layer
    % rock.perm(ucids.unit_cell_ids{33},1) = 700*(milli*darcy);
	% plotToolbar(G, log10(rock.perm(:,1)/(milli*darcy))); axis equal; view([90 0])
    % colormap(turbo)
    % ylim([12500 13500]); zlim([1200 1500])
    
    if strcmp(mesh.type, 'sc2')
        [G, ucids] = updateUpperCentroids(G, ucids, 0);
    end
      
end
% Figures
if fig_perm == 1
    plotRock(G, rock, ucids);
end

%% Helper functions
    function [G, ucids] = updateUpperCentroids(G, ucids, plotFig)
        % NOTE that in case of different grid, cell_ids indices used below
        % (eg 1, 34, 36, 2, 58) must be updated.
        
        % Get center-up reservoir cells
        cell_ids = ucids.unit_cell_ids;
        c_LM2right = cell_ids{1}(~ismember(cell_ids{1}, [cell_ids{34} cell_ids{36}]));
        c_LM2midUp = c_LM2right(all([G.cells.centroids(c_LM2right, 3) < 2325, ...
            G.cells.centroids(c_LM2right, 2) < 27071], 2));
        
        % Find top layer of cells in reservoir  (Ampb neighbors)
        %fc = G.cells.faces(G.cells.facePos(14):G.cells.facePos(14+1)-1,:);
        f2cn = gridCellNo(G);
        fc = G.cells.faces(ismember(f2cn, c_LM2midUp));
        neigh = G.faces.neighbors(fc,:);
        c_Ampb = [cell_ids{2} cell_ids{58}]';
        id_Ampb = [ismember(neigh(:,1), c_Ampb), ...
            ismember(neigh(:,2), c_Ampb)];
        c_LM2top = [neigh(id_Ampb(:,1), 2); neigh(id_Ampb(:,2), 1)];
        id_LM2all = [ismember(neigh(:,1), c_LM2top), ...
            ismember(neigh(:,2), c_LM2top)];
        c_LM2all = [neigh(id_LM2all(:,1), 2); neigh(id_LM2all(:,2), 1)];
        c_LM2all = c_LM2all(ismember(c_LM2all,c_LM2midUp)); % connected to top
        c_LM2all = [c_LM2top; c_LM2all];                    % all
        cn = (G.cells.num - G.layerSize + 1):G.cells.num;
        %c_LM2all_last = cn(ismember(cn, c_LM2all));         % all, last x layer
        c_LM2top_last = cn(ismember(cn, c_LM2top));         % top, last x layer
        
        % Interpolate and assign cell centroids based on y
        zc = interp1(G.cells.centroids(c_LM2top_last, 2), ...
            G.cells.centroids(c_LM2top_last, 3), ...
            G.cells.centroids(c_LM2all, 2));
        gcini = G.cells.centroids(c_LM2all, 3);
        G.cells.centroids(c_LM2all, 3) = zc;
        n = numel(ucids.unit_cell_ids);
        ucids.unit_cell_ids{n+1} = c_LM2top;
        ucids.unit_cell_ids{n+2} = c_LM2all;
        ucids.unitnames{n+1} = 'LM2TopTop';
        ucids.unitnames{n+2} = 'LM2TopAll';
        
        if plotFig
            figure(123)
            plotToolbar(G, ones(G.cells.num,1), 'facecolor','none'); view([90 0]);
            axis equal
            hold on
            plotGrid(G,c_LM2all)
            plot3(repelem(45000,numel(c_LM2all),1), G.cells.centroids(c_LM2all,2), ...
                  gcini,'or')
            plot3(repelem(45000,numel(c_LM2all),1), G.cells.centroids(c_LM2all,2), ...
                G.cells.centroids(c_LM2all,3), '.b')
            
            figure(124)
            plot3(G.cells.centroids(c_LM2all,1), G.cells.centroids(c_LM2all,2), ...
                G.cells.centroids(c_LM2all,3), '.k')
            set(gca,'Zdir','reverse'); view([90 0])
        end
        
    end

    function [rock, ucids, idSealLyrs, idResLyrs] = setLayerPoroPerm(mesh, rock, ucids, opt)
        unit_cells = ucids.unit_cell_ids;
        if opt.sc == 1
            idSealLyrs = 2:22;  
            idResLyrs  = [];
        elseif opt.sc == 2 && strcmp(opt.studyCase, 'sc2_baseCase') 
            idSealLyrs = 2:2:22; 
            idResLyrs  = 3:2:21;
        elseif opt.sc == 2 && strcmp(opt.studyCase, 'sc2_predict_theta30_coarse') ...
               || opt.sc == 2 && strcmp(opt.studyCase, 'sc2_predict_theta30') ...
               || opt.sc == 2 && strcmp(opt.studyCase, 'sc2_predict_theta30_coarse_nohyst') ...
               || opt.sc == 2 && strcmp(opt.studyCase, 'sc2_predict_theta30_nohyst')
            idSealLyrs = 2:2:22; 
            idResLyrs  = 3:2:21;
        else
            error('only scenario 1 (continuous seal) and 2 (disc.) are possible')
        end
        idu = [1:(1+mesh.ampLyr+2) 58];
        idfu = idu(1:end-1);
        idamp = 58; idmio = 23; idyou = 24;
%         switch mesh.type
%             case 'sc2'
%                 if mesh.ampLyr == 21 && ucids.upper == 1
%                     idu = [1:(1+mesh.ampLyr+2) 58];
%                     idfu = idu(1:end-1);
%                     idamp = 58; idmio = 23; idyou = 24;
%                     idSealLyrs = 2:2:22; idResLyrs = 3:2:21;
%                 end
%                 
%             otherwise
%                 if strcmp(mesh.type, 'coarse') || strcmp(mesh.type, 'ref300') || ...
%                         strcmp(mesh.type, 'ref200') || strcmp(mesh.type, 'ref100') || ...
%                         strcmp(mesh.type, 'ref50') || strcmp(mesh.type, 'ref40')
%                     idu  = 1:7;
%                     idfu = 4:7;
%                     
%                     % Anahuac and older (Oligocene and older; all considered Anahuac here)
%                     k_h = 0.001;                                % Lu et al (2015)
%                     rock.poro(unit_cells{1}, 1) = 0.08;         % "
%                     rock.perm(unit_cells{1}, [1, 4]) = k_h*milli*darcy;
%                     rock.perm(unit_cells{1}, 6) = 0.2*k_h*milli*darcy;
%                     % LM1 reservoir
%                     k_h = 100;                                  % Wallace (2013); Wallace et al. (2017, GOM atlas);
%                     poro = log(k_h/.7385)/20.011;               % "
%                     rock.poro(unit_cells{2}, 1) = poro;
%                     rock.perm(unit_cells{2}, [1, 4]) = k_h*milli*darcy;
%                     rock.perm(unit_cells{2}, 6) = 0.2*k_h*milli*darcy;
%                     % Marg A
%                     k_h = 0.003;                                % Lu et al. (2017, GOM atlas)
%                     rock.poro(unit_cells{3}, 1) = 0.1;          % "
%                     rock.perm(unit_cells{3}, [1, 4]) = k_h*milli*darcy;
%                     rock.perm(unit_cells{3}, 6) = (1/3)*k_h*milli*darcy;
%                     
%                 elseif strcmp(mesh.type, 'fine')
%                     % later
%                     
%                 elseif ucids.upper == 1    % only upper section
%                     idu  = 1:4;
%                     idfu = idu;
%                 elseif ucids.upper == 2    % only resLM2 and amp
%                     idu = 1:2;
%                     idfu = idu;
%                 end
%         end
        
        % unit ids
        ucids.id.units         = idu; 
        ucids.id.faulted_units = idfu;                     % unit ids
        ucids.id.fault         = numel(unit_cells);        % fault ids
        if strcmp(mesh.type, 'sc2') && mesh.ampLyr == 21
           ucids.fault.juxt = [idfu(1), idfu(1); idfu(1), idfu(2);  idfu(1), idfu(3); idfu(1), idfu(4); ...
                               idfu(1), idfu(5); idfu(2), idfu(6); idfu(3), idfu(7); idfu(4), idfu(8); ...
                               idfu(5), idfu(9); idfu(6), idfu(10); idfu(7), idfu(11); idfu(8), idfu(12); ...
                               idfu(9), idfu(13); idfu(10), idfu(14); idfu(11), idfu(15); idfu(12), idfu(16); ... 
                               idfu(13), idfu(17); idfu(14), idfu(18); idfu(15), idfu(19); idfu(16), idfu(20);
                               idfu(17), idfu(21); idfu(18), idfu(22); idfu(19), idfu(23); idfu(20), idfu(23);
                               idfu(21), idfu(23); idfu(22), idfu(23); idfu(23), idfu(23); idfu(24), idfu(24)];
                               
%         elseif isfield(ucids, 'upper') && ucids.upper == 2
%             ucids.fault.juxt = [idfu(1), idfu(1); idfu(1), idfu(2); ...
%                                 idfu(2), idfu(2)];
%         else
%             ucids.fault.juxt = [idfu(1), idfu(1); idfu(1), idfu(2); ...
%                                 idfu(2), idfu(2); idfu(2), idfu(3); ...
%                                 idfu(3), idfu(3); idfu(4), idfu(4)];
        end
        
        % Vcl (see fig 4.12 in Trevino & Meckel, Atlas GoM, 2017)
        grlog = sortrows(load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/data_files/Well1.dat'), 2);
        grlog(:,2) = grlog(:,2)*ft;
        grmM = [prctile(grlog(:,1), 10); prctile(grlog(:,1), 90)];         
        Vsh = (grlog(:,1)-grmM(1)) ./ (grmM(2) - grmM(1));
        Vsh(Vsh>1) = 1;     Vsh(Vsh<0) = 0;
        lims = [7079.5 8049.3 10196 10622]*ft;
        uids = [grlog(:,2)<=lims(1), all([grlog(:,2)>lims(1), grlog(:,2)<=lims(2)], 2), ...
                all([grlog(:,2)>lims(2), grlog(:,2)<=lims(3)], 2), ...
                all([grlog(:,2)>lims(3), grlog(:,2)<=lims(4)], 2), ...
                grlog(:,2)>lims(4)];
       ucids.Vcl = [];
       VshToVcl = 0.65;                                                     % transform Vsh to Vcl
       if strcmp(mesh.type, 'sc2') && mesh.ampLyr==21
           switch opt.petr{3}
               case 'sandLogAvg_clayLogAvg'
                   ucids.Vcl(1:58) = nan;
                   ucids.Vcl(1) = VshToVcl*mean(Vsh(uids(:,3)));           % res LM2
                   ucids.Vcl(idSealLyrs) = VshToVcl*mean(Vsh(uids(:,2)));  % Amph B layering seal
                   ucids.Vcl(idResLyrs) = VshToVcl*mean(Vsh(uids(:,2)));   % Amph b layering reserv
                   ucids.Vcl(58) = mean([ucids.Vcl(idSealLyrs) ucids.Vcl(idResLyrs)]);   % Amph b rest
                   ucids.Vcl(23) = VshToVcl*mean(Vsh(uids(:,1)));          % MM-UM
                   ucids.Vcl(24) = 0;                                      % Younger
                   
               case 'sandLogMin_clayLog70'
                   ucids.Vcl(1:58) = nan;
                   ucids.Vcl(1) = VshToVcl*mean(Vsh(uids(:,3)));                        % res LM2, consider average since no resolution
                   ucids.Vcl(idSealLyrs) = VshToVcl*prctile(Vsh(uids(:,2)), 70);        % Amph B layering seal
                   ucids.Vcl(idResLyrs) = VshToVcl*min(Vsh(uids(:,2)));                 % Amph b layering reserv
                   ucids.Vcl(58) = mean([ucids.Vcl(idSealLyrs) ucids.Vcl(idResLyrs)]);  % Amph b rest
                   ucids.Vcl(23) = VshToVcl*mean(Vsh(uids(:,1)));                       % MM-UM
                   ucids.Vcl(24) = 0;                                                   % Younger
                   
               case 'sand0_clayLogAvg'
                   ucids.Vcl(1:58) = nan;
                   ucids.Vcl(1) = 0;                                       % res LM2
                   ucids.Vcl(idSealLyrs) = VshToVcl*mean(Vsh(uids(:,2)));  % Amph B layering seal
                   ucids.Vcl(idResLyrs) = 0;                               % Amph b layering reserv
                   ucids.Vcl(58) = mean([ucids.Vcl(idSealLyrs) ucids.Vcl(idResLyrs)]);         % Amph b rest
                   ucids.Vcl(23) = 0;                                      % MM-UM
                   ucids.Vcl(24) = 0;                                      % Younger
           end
       end
%                 elseif strcmp(mesh.type, 'coarse') || strcmp(mesh.type, 'ref300') || ...
%                    strcmp(mesh.type, 'ref200') || strcmp(mesh.type, 'ref100') || ...
%                    strcmp(mesh.type, 'ref50') || strcmp(mesh.type, 'ref40')
%                     ucids.Vcl(1) = 0.8;                                     % Anahuac and older
%                     ucids.Vcl(2) = 0;                                       % LM1 reservoir
%                     ucids.Vcl(3) = mean(Vsh(uids(:,4)));                    % Marg A
%                     ucids.Vcl(4) = 0;                                       % res LM2
%                     ucids.Vcl(5) = mean(Vsh(uids(:,2)));                    % Amph B
%                     ucids.Vcl(6) = 0;                                       % MM-UM
%                     ucids.Vcl(7) = 0;                                       % Younger
%                 elseif ucids.upper == 1
%                     ucids.Vcl(1) = 0;                                       % res LM2
%                     ucids.Vcl(2) = mean(Vsh(uids(:,2)));                    % Amph B
%                     ucids.Vcl(3) = 0;                                       % MM-UM
%                     ucids.Vcl(4) = 0;                                       % Younger
%                 elseif ucids.upper == 2
%                     ucids.Vcl(1) = 0;                                       % res LM2
%                     ucids.Vcl(2) = mean(Vsh(uids(:,2)));                    % Amph B
%                end
%            case 'sandLogAvg_clayLogAvg'
%                if strcmp(mesh.type, 'coarseUpper')
%                     ucids.Vcl(1) = mean(Vsh(uids(:,3)));                    % res LM2
%                     ucids.Vcl(2) = mean(Vsh(uids(:,2)));                    % Amph B
%                     ucids.Vcl(3) = mean(Vsh(uids(:,1)));                    % MM-UM
%                     ucids.Vcl(4) = 0;                                       % Younger
%                end
%            case 'sand1_clayLogAvg'
%                if strcmp(mesh.type, 'coarseUpper')
%                     ucids.Vcl(1) = 0.1;                                     % res LM2
%                     ucids.Vcl(2) = mean(Vsh(uids(:,2)));                    % Amph B
%                     ucids.Vcl(3) = 0.1;                                     % MM-UM
%                     ucids.Vcl(4) = 0;                                       % Younger
%                end
%       end
        
       % Amph B
       k_ha = 0.005;                                      % Lu et al. (2017, GOM atlas)
       %if strcmp(mesh.type, 'sc2')
       % Sealing sublayers
       poros = 0.15;
       rock.poro([unit_cells{idSealLyrs}], 1) = poros;        % "
       rock.perm([unit_cells{idSealLyrs}], [1, 4]) = k_ha*milli*darcy;
       rock.perm([unit_cells{idSealLyrs}], 6) = (1/5)*k_ha*milli*darcy;
       
       % Permeable sublayers
       k_hap = 175; % mD
       poro = log(k_hap/.7385)/20.011;
       rock.poro([unit_cells{idResLyrs}], 1) = poro;
       rock.perm([unit_cells{idResLyrs}], [1, 4]) = k_hap*milli*darcy;
       rock.perm([unit_cells{idResLyrs}], 6) = (1/3)*k_hap*milli*darcy;
       
       % Amph b rest (average)
       if opt.sc == 1
           rock.poro(unit_cells{idamp}, 1) = poros;        % "
           rock.perm(unit_cells{idamp}, [1, 4]) = k_ha*milli*darcy;
           rock.perm(unit_cells{idamp}, 6) = (1/5)*k_ha*milli*darcy;
       elseif opt.sc == 2
           rock.poro(unit_cells{idamp}, 1) = mean([poros poro]);        % "
           rock.perm(unit_cells{idamp}, [1, 4]) = mean([k_ha k_hap])*milli*darcy;
           rock.perm(unit_cells{idamp}, 6) = harmmean([(1/5)*k_ha (1/3)*k_hap])*milli*darcy;
       end
           
           % MM-UM
           k_hm = 200;                                       % Wallace (2013); Wallace et al. (2017, GOM atlas);
           poro = log(k_hm/.7385)/20.011;                    % "
           rock.poro(unit_cells{idmio}, 1) = poro;
           rock.perm(unit_cells{idmio}, [1, 4]) = k_hm*milli*darcy;
           rock.perm(unit_cells{idmio}, 6) = (1/3)*k_hm*milli*darcy;
           
           % Younger
           k_hy = 500;                                     % Guess (not that important here)
           rock.poro(unit_cells{idyou}, 1) = 0.35;         %
           rock.perm(unit_cells{idyou}, [1, 4]) = k_hy*milli*darcy;
           rock.perm(unit_cells{idyou}, 6) = k_hy*milli*darcy;
           
%        else
%            rock.poro(unit_cells{idfu(2)}, 1) = 0.15;        % "
%            rock.perm(unit_cells{idfu(2)}, [1, 4]) = k_ha*milli*darcy;
%            rock.perm(unit_cells{idfu(2)}, 6) = (1/3)*k_ha*milli*darcy;
%            if ~isfield(ucids, 'upper') || ucids.upper ~= 2
%                % MM-UM
%                k_h = 200;                                       % Wallace (2013); Wallace et al. (2017, GOM atlas);
%                poro = log(k_h/.7385)/20.011;                    % "
%                rock.poro(unit_cells{idfu(3)}, 1) = poro;
%                rock.perm(unit_cells{idfu(3)}, [1, 4]) = k_h*milli*darcy;
%                rock.perm(unit_cells{idfu(3)}, 6) = (1/3)*k_h*milli*darcy;
%                % Younger
%                k_h = 500;                                       % Guess (not that important here)
%                rock.poro(unit_cells{idfu(4)}, 1) = 0.35;        %
%                rock.perm(unit_cells{idfu(4)}, [1, 4]) = k_h*milli*darcy;
%                rock.perm(unit_cells{idfu(4)}, 6) = k_h*milli*darcy;
%            end
%        end
        
        % res LM2
        if strcmp(mesh.type, 'sc2')
            cells_resLM2 = unit_cells{1}; 
        else
            cells_resLM2 = unit_cells{idfu(1)};
        end
        k_h = 150; % mD
        k_v = (1/5)*k_h;
        poro = log(k_h/.7385)/20.011;
        %S_wr = log(k_h/4025.1)/-10.64;
        switch opt.petr{1}
            case 'isotr-homog' % isotropic and homogeneous for all cells in layer.
                rock.poro(cells_resLM2, 1) = poro;
            case 'isotr-rand'  % isotropic within each cell, random for each cell.
                [poromin,poromax] = deal(poro-0.1,poro+0.1);
                porovals = (poromax-poromin).*rand(numel(cells_resLM2),1) + poromin;
                rock.poro(cells_resLM2, 1) = porovals;
            otherwise
                error('porosity case not supported.')
        end
        
        switch opt.petr{2}
            case 'isotr-homog' % isotropic and homogeneous for all cells in layer.
                rock.perm(cells_resLM2, [1,4,6]) = k_h*milli*darcy;
            case 'isotr-rand'  % isotropic within each cell, random for each cell.
                [kmin,kmax] = deal(k_h-100,k_h+100);
                kvals = (kmax-kmin).*rand(numel(cells_resLM2),1) + kmin;
                rock.perm(cells_resLM2, [1, 4]) = kvals*milli*darcy;
                rock.perm(cells_resLM2, 6) = kvals*milli*darcy;
            case 'anisotropic' % kx = 5kz (Wallace, 2013; Wallace et al., 2017)
                rock.perm(cells_resLM2, [1, 4]) = k_h*milli*darcy;
                rock.perm(cells_resLM2, 6) = k_v*milli*darcy;
            otherwise
                error('permeability case not supported.')
        end
    end

    function [rock, fault] = setFaultPoroPerm(mesh, rock, ucids, opt, G, tlim, idSealLyrs, idResLyrs, fig_perm)
        unit_cells = ucids.unit_cell_ids;
        switch opt.fk_pred
            case 'test'
                % Assign basal fault perm and poro
                rock.poro(unit_cells{end}, 1) = 0.1;
                rock.perm(unit_cells{end}, [1, 4]) = 1e-3*milli*darcy;
                rock.perm(unit_cells{end}, 6) = 1e-2*milli*darcy;

            case 'predict'
                % GRID DEPENDENT!
                fault.fcells = [unit_cells{26:31}]'; % not sorted

                % Poro: UPTDATE with predict porosity later
                fault.poro = 0.15*ones(numel(fault.fcells), 1);
                
                % Perm: from PREDICT parametrized distribution
                % fault is 26 to 31 (famp1 to famp6), the rest of fault
                % cells we will use the SGR.
                fault.k = zeros(numel(fault.fcells), 3);
                ups_pth = fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/upscaling/upscaled_data/');
                fnames = {'famp1_kfit', 'famp2_kfit', 'famp3_kfit', ...
                          'famp4_kfit', 'famp5_kfit', 'famp6_kfit'};
                xctr = unique(G.cells.centroids(:,1));
                ncy = numel(mesh.thick);
                assert(numel(xctr) == ncy);
                idf = 26;
                for n=1:numel(fnames)
                    kdat = load([ups_pth fnames{n} '.mat']);
                    kvals =  [random(kdat.pdx, [ncy,1]), random(kdat.pdy, [ncy,1]), ...
                              random(kdat.pdz, [ncy,1])]; % log10_mD
                    kvals = 10.^(kvals)*(milli*darcy);
                    for j = 1:ncy
                        cint = unit_cells{idf};
                        id_cint = G.cells.centroids(cint,1)==xctr(j);
                        cid = cint(id_cint);
                        fault.k(ismember(fault.fcells, cid), :) = ...
                                        repmat(kvals(j,:), numel(cid), 1);
                    end
                    idf = idf + 1;
                end

                % Transform perm and assign to rock
                rock.poro(fault.fcells, 1) = fault.poro;
                switch mesh.type
                    case 'sc2'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/sc2/throwDiv/fnodcoord.mat');
                end
                [rock, fault] = transformPerm(G, rock, fault, 'rotateAboutX');
                % plotCellData(G,log10(rock.perm(:,end)/(milli*darcy)), fault.fcells, 'edgecolor', 'none')

                % Add LM2, Miocene, Young from SGR run
                if strcmp(opt.studyCase, 'sc2_predict_theta30_nohyst') ...
                   || strcmp(opt.studyCase, 'sc2_predict_theta30')
                    fSGR = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/data_files/faultPermSGR_Gupper63lyr.mat');
                else
                    fSGR = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/data_files/faultPermSGR_Gupper29lyr.mat');
                end
                rock.perm(unit_cells{25},:) = fSGR.rock.perm(unit_cells{25},:); %fLM2
                rock.perm(unit_cells{32},:) = fSGR.rock.perm(unit_cells{32},:); %fMMUM
                rock.perm(unit_cells{33},:) = fSGR.rock.perm(unit_cells{33},:); %fYoung
                rock.poro(unit_cells{25}) = fSGR.rock.poro(unit_cells{25});
                rock.poro(unit_cells{32}) = fSGR.rock.poro(unit_cells{32});
                rock.poro(unit_cells{33}) = fSGR.rock.poro(unit_cells{33});

                % plotCellData(G,log10(rock.perm(:,1)/(milli*darcy)), unit_cells{end}, 'edgecolor', 'none')
                % colormap(turbo); colorbar; view([60 -10])

            case 'spe02'
                predictor      = 'SGR';
                ucids.juxt     = 0;
                G              = cellUnitTags(G, ucids, 0);
                optSmooth      = false;
                optContacts    = true;
                [kPred, throw, ncont] = faultPermPred(G, ucids, predictor, optSmooth, optContacts, tlim);
                G.cells.zmax   = G.cells.centroids(:, 3);
                fault.zf       = throw;                     
                fault.zf(fault.zf > 50) = 50;                              % 50m assumed same everywhere.
                fault.throw    = throw;
                fault.fcells   = unit_cells{end}';
                if isfield(opt, 'SGR') && strcmp(opt.SGR, 'manual')
                    ids = ismember(fault.fcells, [unit_cells{27:30}]);     % for 25m subdivisions in caprock
                    maxV = mean([mean(ucids.Vcl(idSealLyrs)) mean(ucids.Vcl(idResLyrs))]);
                    kPred.SGR(ids) = maxV;
                    kPred.SGR(kPred.SGR > maxV) = maxV;
                end
                fault.kPred    = kPred;
                fault.k        = predToPerm(G, rock, fault, 'name', 'spe02', ...
                                            'kAlong', true, 'kAniso', 10, ...
                                            'cutoff', '200', 'smooth', false);
                [sPoro, cPoro] = endMemberPoro(G, rock, ucids, 'rhoBrine', 1.0689e+03, ...
                                               'sandOpt', opt.fSandPoro);
                fault.poro     = relatePoroPerm(G, rock, fault, 'out', 'poro', ...
                                                'method', 'idealPack', ...
                                                'endSand', sPoro, 'endClay', cPoro);
                assignin('base', 'ncont', ncont);
                
                % Assign properties to rock
                rock.poro(fault.fcells, 1) = fault.poro;
                switch mesh.type
                    case 'coarse'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/fnodcoord.mat'); % updated
                    case 'coarseUpper'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/fnodcoord.mat');
                    case 'fine'
                        % later
                    case 'ref300'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref300/fnodcoord.mat');
                    case 'ref200'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref200/fnodcoord.mat');
                    case 'ref100'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref100/fnodcoord.mat');
                    case 'ref50'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref50/fnodcoord.mat');
                    case 'ref25'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref25/fnodcoord.mat');
                    case 'ref12.5'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref12.5/fnodcoord.mat');
                    case 'sc2'
                        fault.fn = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/sc2/throwDiv/fnodcoord.mat');
                end
                [rock, fault] = transformPerm(G, rock, fault, 'rotateAboutX');    

                % save if needed for predict
                %save('faultPermSGR_Gupper63lyr.mat','fault','rock','ncont')
                   
            case 'smearsUpscaled'                
                
        end
        
        if fig_perm == 1 % Plot fault
            if strcmp(opt.fk_pred,'spe02')
                %cmap = cmap_agu();
                cmap = copper;
                colormap(cmap)
                figure(78)
                kPlot.kPred = zeros(G.cells.num, 1);  kPlot.kFault = zeros(G.cells.num, 1);
                kPlot.porof = zeros(G.cells.num, 1); kPlot.throw = zeros(G.cells.num, 1);
                kPlot.kPred(unit_cells{end}) = kPred.SGR;
                kPlot.kFault(unit_cells{end}) = log10(fault.k.vals(:, 1)/(milli*darcy));
                kPlot.porof(unit_cells{end}) = fault.poro;
                kPlot.throw(unit_cells{end}) = fault.throw;
                cc = max(G.faces.centroids(:,1)/2); limm = 2000;
                idplot1 = G.cells.centroids(unit_cells{end}, 1)>= cc-limm;
                idplot2 = G.cells.centroids(unit_cells{end}, 1) <= cc+limm;
                idplot = all([idplot1, idplot2], 2);
                plotToolbar(G, kPlot, unit_cells{end}(idplot), 'FaceAlpha', 0.7)
                hold on
                col = ['k','w'];
                for n = [1 numel(ncont)-1 numel(ncont)]
                    idf1 = all([G.nodes.coords(ncont{n}(:,1), 1) >= cc-1.1*limm, ...
                        G.nodes.coords(ncont{n}(:,1), 1) <= cc+1.1*limm], 2);
                    idf2 = all([G.nodes.coords(ncont{n}(:,2), 1) >= cc-1.1*limm, ...
                        G.nodes.coords(ncont{n}(:,2), 1) <= cc+1.1*limm], 2);
                    x = G.nodes.coords(ncont{n}(idf1,1), 1); x2 = G.nodes.coords(ncont{n}(idf2,2), 1);
                    y = G.nodes.coords(ncont{n}(idf1,1), 2); y2 = G.nodes.coords(ncont{n}(idf2,2), 2);
                    z = G.nodes.coords(ncont{n}(idf1,1), 3); z2 = G.nodes.coords(ncont{n}(idf2,2), 3);
                    plot3(x, y, z, col(1), 'lineStyle', '-', 'lineWidth', 2)
                    plot3(x2, y2, z2, col(2), 'lineStyle', '-', 'lineWidth', 2)
                end
                axis equal off
                view([180, 0])
                colormap(copper)
            end
        end
        
    end

    function plotRock(G, rock, ucids)
        ucids = ucids.unit_cell_ids;
        fc = ucids{61};
        h = figure(21);
        plotToolbar(G, log10(rock.perm(:,1)/(milli*darcy)), fc, ...
                    'EdgeColor', [0.5, 0.5, 0.5])
        c = colorbar;
        view([0 0]) %view([20 -20])
        colormap(copper)
        c.Label.String = '$\log_{10}k_\mathrm{xx}$ [mD])';
        c.Label.Interpreter = 'latex';
        c.Label.FontSize = 16;
        %title('log10(k_{xx} [mD])')
        set(gca, 'fontSize', 12)
        xv = [15000 30000];
        zv = [1200 2200];
        xlim(xv); zlim(zv);
        % --- Adjust caxis to plotting area 
        idz = all([G.cells.centroids(fc,end) < max(zv), ...
                   G.cells.centroids(fc,end) > min(zv)], 2);
        idx = all([G.cells.centroids(fc,1) < max(xv), ...
                   G.cells.centroids(fc,1) > min(xv)], 2);  
        idf = all([idz, idx], 2);
        min_logk = floor(min(log10(rock.perm(fc(idf),1)/(milli*darcy))));
        max_logk = ceil(max(log10(rock.perm(fc(idf),1)/(milli*darcy))));
        caxis([min_logk max_logk]);
        % ---
        axis off
        xlabel('$x$ [m]', 'interpreter', 'latex', 'fontsize', 12)
        ylabel('$y$ [m]', 'interpreter', 'latex', 'fontsize', 12)
        zlabel('$z$ [m]', 'interpreter', 'latex', 'fontsize', 12)
        set(h, 'position', [50, 50, 1000, 500])
        
        figure(22)
        plotCellData(G, rock.poro, [ucids{end-6} ucids{end-5} ...
                     ucids{end-4} ucids{end-3}], ...
                     'EdgeColor', [0.5, 0.5, 0.5])
        axis equal; c = colorbar; caxis([0 1]);
        set(gca,'Xdir','reverse'); axis equal off; view([90 0])
        
        figure(23)
        %plotCellData(G, rock.perm(:,1)/(milli*darcy), ucids{4}, 'EdgeColor', [0.6 0.6 0.6])
        %plotCellData(G, rock.perm(:,1)/(milli*darcy), ucids{5}, 'EdgeColor', 'w')
        plotCellData(G, log10(rock.perm(:,1)/(milli*darcy)), 'EdgeColor', 'none')
        axis equal; colormap(copper); colorbar
        ylim([10000, 20000]); set(gca,'Xdir','reverse'); axis equal off; view([90 0])
        
        figure(24)
        plotCellData(G, rock.poro, [ucids{end-10} ucids{end-9}], 'EdgeColor', 'none')
        axis equal; c = colorbar; caxis([0 1]);
        ylim([10000, 16000]); set(gca,'Xdir','reverse'); axis equal off; view([90 0])
        
        figure(25)
        cmap = [236, 214, 56; 102, 77, 0]./255;
        colormap(cmap);
        idx = ones(G.cells.num,1); 
        idx([ucids{1} ucids{end-7} ucids{end-6} ucids{end-5} ucids{end-4} ucids{end-3}]) = 2;
        plotToolbar(G, idx), set(gca,'Xdir','reverse');
        view(55, 25); camproj perspective; axis equal tight 
        caxis([1 2])
        
        figure(26)
        cmap = [236, 214, 56; 226, 204, 46; 102, 77, 0; 125 125 125]./255;
        colormap(cmap);
        idx = ones(G.cells.num, 1); 
        idx([ucids{end-2}]) = 2;
        idx([ucids{end-5} ucids{end-3}]) = 3;
        idx([ucids{end-4}]) = 4;
        cc = max(G.faces.centroids(:,1)/2);
        idplot1 = G.cells.centroids(ucids{end}, 1)>= cc-30;
        idplot2 = G.cells.centroids(ucids{end}, 1) <= cc+30;
        idplot = all([idplot1, idplot2], 2);
        plotToolbar(G, idx, ucids{end}(idplot)), set(gca,'Xdir','reverse');
        view(55, 25); camproj perspective; axis equal off 
        caxis([1 4])
        
        %figure(26)
        %plotCellData(G, rock.poro, 1:G.layerSize)
        %c = colorbar; caxis([0 1]);
        %view([70 30]); axis off; axis equal
    end

%% Preliminary step: Obtain cell ids of each unit
if obtained.ucids == 0
    switch mesh.type
        case 'coarse'
            dirc = 'mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/cell_centroids_';
            unitnames = {'old', 'resLM1', 'marg', 'resLM2', 'amp', 'MMUM', ...
                         'you', 'fLM2', 'fLM2_amp', 'famp', 'famp_MMUM', ...
                         'fMMUM', 'fyou'}';
            centroidsL = dlmread('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/resLM2LeftBlock/cell_centroids_resLM2Left.dat'); 
        case 'coarseUpper'
            dirc = 'mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/cell_centroids_';
            unitnames = {'resLM2', 'amp', 'MMUM', 'you', 'fLM2', 'fLM2_amp', 'famp', 'famp_MMUM', ...
                         'fMMUM', 'fyou'}';
            centroidsL = dlmread('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/resLM2LeftBlock/cell_centroids_resLM2Left.dat');
        case 'fine'
            % later
            
        case 'ref300'
            dirc = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref300/cell_centroids_';
            unitnames = {'old', 'resLM1', 'marg', 'resLM2', 'amp', 'MMUM', 'you', 'fLM2', 'fLM2_amp', 'famp', 'famp_MMUM', ...
                         'fMMUM', 'fyou'}';
            centroidsRef = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref300/cell_centroids_resLM2ref.dat');
            centroidsL = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref300/cell_centroids_resLM2Left.dat');
            elemOnTop  = transpose(sort(dlmread('mrst-dev/mrst-ls/ls-proj/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref300/elemOnTopAreaInWedge_2D.dat')));            
        case 'ref200'
            dirc = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref200/cell_centroids_';
            unitnames = {'old', 'resLM1', 'marg', 'resLM2', 'amp', 'MMUM', 'you', 'fLM2', 'fLM2_amp', 'famp', 'famp_MMUM', ...
                         'fMMUM', 'fyou'}';
            centroidsRef = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref200/cell_centroids_resLM2ref.dat');
            centroidsL = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref200/cell_centroids_resLM2Left.dat');
            elemOnTop  = transpose(sort(dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref200/elemOnTopAreaInWedge_2D.dat')));      
        case 'ref100'
            dirc = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref100/cell_centroids_';
            unitnames = {'old', 'resLM1', 'marg', 'resLM2', 'amp', 'MMUM', 'you', 'fLM2', 'fLM2_amp', 'famp', 'famp_MMUM', ...
                         'fMMUM', 'fyou'}';
            centroidsRef = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref100/cell_centroids_resLM2ref.dat');
            centroidsL = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref100/cell_centroids_resLM2Left.dat');
            elemOnTop  = transpose(sort(dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref100/elemOnTopAreaInWedge_2D.dat')));
        case 'ref50'
            dirc = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref50/cell_centroids_';
                 unitnames = {'old', 'resLM1', 'marg', 'resLM2', 'amp', 'MMUM', 'you', 'fLM2', 'fLM2_amp', 'famp', 'famp_MMUM', ...
                              'fMMUM', 'fyou'}';
            centroidsRef = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref50/cell_centroids_resLM2ref.dat');
            centroidsL = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref50/cell_centroids_resLM2Left.dat');
            elemOnTop  = transpose(sort(dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref50/elemOnTopAreaInWedge_2D.dat')));
        case 'ref25'
            dirc = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref25/cell_centroids_';
            if mesh.reduce == 1
                unitnames = {'resLM2', 'amp', 'MMUM', 'you', 'fLM2', 'fLM2_amp', 'famp', 'famp_MMUM', ...
                    'fMMUM', 'fyou'}';
            else
                unitnames = {'old', 'resLM1', 'marg', 'resLM2', 'amp', 'MMUM', 'you', 'fLM2', 'fLM2_amp', 'famp', 'famp_MMUM', ...
                    'fMMUM', 'fyou'}';
            end
            centroidsRef = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref25/cell_centroids_resLM2ref.dat');
            centroidsL = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref25/cell_centroids_resLM2Left.dat');
            elemOnTop  = transpose(sort(dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref25/elemOnTopAreaInWedge_2D_upper.dat')));
       case 'ref12.5'
            dirc = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref12.5/cell_centroids_';
            if mesh.reduce == 2
                unitnames = {'resLM2', 'amp', 'fLM2', 'fLM2_amp', 'famp'}';
            else
                unitnames = {'old', 'resLM1', 'marg', 'resLM2', 'amp', 'MMUM', 'you', 'fLM2', 'fLM2_amp', 'famp', 'famp_MMUM', ...
                    'fMMUM', 'fyou'}';
            end
            centroidsRef = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref12.5/cell_centroids_resLM2ref.dat');
            centroidsL = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/gridRef/ref12.5/cell_centroids_resLM2Left.dat');
            elemOnTop  = transpose(sort(dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref12.5/elemOnTopAreaInWedge_2D_upper.dat')));
        case 'sc21'
            dirc = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/scenario2/initial/cell_centroids_';
            if mesh.reduce == 1
                unitnames = {'res_LM2',  'AmphB', 'MMUM', 'Younger', ...
                             'f_LM2', 'f_amp', 'f_MMUM', 'f_young', 'Ar1', 'Ar2', 'Ar3', ...
                             'Ar4', 'Ar5', 'Ar6', 'Ar7', 'Ar8', 'Ar9', 'Ar10', 'Ar11', 'Ar12', ...
                             'Ar13', 'Ar14', 'Ar15', 'Ar16', 'Ar17', 'Ar18', 'Ar19', 'Ar20', 'Ar21', ...
                             'Al1', 'Al2', 'Al3', 'Al4', 'Al5', 'Al6', 'Al7', 'Al8', 'Al9', 'Al10', ...
                             'Al11', 'Al12', 'Al13', 'Al14', 'Al15', 'Al16', 'Al17', 'Al18', 'Al19', ...
                             'Al20', 'Al21', 'res_LM2_left', 'res_LM2_ref', 'res_LM2_refl', ...
                             'AmphB_left', 'MMUM_left'}';
            else
                unitnames = {'res_LM2', 'AmphB', 'MMUM', 'Younger', ...
                             'f_LM2', 'f_amp', 'f_MMUM', 'f_young', 'Ar1', 'Ar2', 'Ar3', ...
                             'Ar4', 'Ar5', 'Ar6', 'Ar7', 'Ar8', 'Ar9', 'Ar10', 'Ar11', 'Ar12', ...
                             'Ar13', 'Ar14', 'Ar15', 'Ar16', 'Ar17', 'Ar18', 'Ar19', 'Ar20', 'Ar21', ...
                             'Al1', 'Al2', 'Al3', 'Al4', 'Al5', 'Al6', 'Al7', 'Al8', 'Al9', 'Al10', ...
                             'Al11', 'Al12', 'Al13', 'Al14', 'Al15', 'Al16', 'Al17', 'Al18', 'Al19', ...
                             'Al20', 'Al21',  'Older', 'res_LM1', 'MargA', 'res_LM2_left', 'res_LM2_ref', 'res_LM2_refl', ...
                             'AmphB_left', 'MMUM_left'}';
            end
            centroidsRef = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/scenario2/initial/cell_centroids_res_LM2_ref.dat');
            centroidsLref = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/scenario2/initial/cell_centroids_res_LM2_refl.dat');
            centroidsL = dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/scenario2/initial/cell_centroids_res_LM2_left.dat');
            elemOnTop  = transpose(sort(dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref12.5/elemOnTopAreaInWedge_2D_upper.dat')));  
            
            case 'sc2'
            dirc = 'mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/output_treltogprs/extr_tri/2D/scenario2/throwDiv/cell_centroids_';
            if mesh.reduce == 1
                unitnames = {'res_LM2',  'Ar1', 'Ar2', 'Ar3', ...
                             'Ar4', 'Ar5', 'Ar6', 'Ar7', 'Ar8', 'Ar9', 'Ar10', 'Ar11', 'Ar12', ...
                             'Ar13', 'Ar14', 'Ar15', 'Ar16', 'Ar17', 'Ar18', 'Ar19', 'Ar20', 'Ar21', ...
                             'MMUM', 'Younger', 'f_LM2', 'f_amp1', 'f_amp2', 'f_amp3', 'f_amp4', 'f_amp5', 'f_amp6', ...
                             'f_MMUM', 'f_young', 'res_LM2_left', 'res_LM2_ref', 'res_LM2_refl', ...
                             'Al1', 'Al2', 'Al3', 'Al4', 'Al5', 'Al6', 'Al7', 'Al8', 'Al9', 'Al10', ...
                             'Al11', 'Al12', 'Al13', 'Al14', 'Al15', 'Al16', 'Al17', 'Al18', 'Al19', ...
                             'Al20', 'Al21', 'AmphB', 'AmphB_left', 'MMUM_left'}';
            else
                unitnames = {'res_LM2', 'AmphB', 'MMUM', 'Younger', ...
                             'f_LM2', 'f_amp1', 'f_amp2', 'f_amp3', 'f_amp4', 'f_amp5', 'f_amp6', ...
                             'f_MMUM', 'f_young', 'Ar1', 'Ar2', 'Ar3', ...
                             'Ar4', 'Ar5', 'Ar6', 'Ar7', 'Ar8', 'Ar9', 'Ar10', 'Ar11', 'Ar12', ...
                             'Ar13', 'Ar14', 'Ar15', 'Ar16', 'Ar17', 'Ar18', 'Ar19', 'Ar20', 'Ar21', ...
                             'Al1', 'Al2', 'Al3', 'Al4', 'Al5', 'Al6', 'Al7', 'Al8', 'Al9', 'Al10', ...
                             'Al11', 'Al12', 'Al13', 'Al14', 'Al15', 'Al16', 'Al17', 'Al18', 'Al19', ...
                             'Al20', 'Al21', 'Older', 'res_LM1', 'MargA', 'res_LM2_left', 'res_LM2_ref', 'res_LM2_refl', ...
                             'AmphB_left', 'MMUM_left'}';
            end
            elemOnTop  = transpose(sort(dlmread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/sc2/throwDiv/elemOnTopAreaInWedge_2D_upper.dat')));
    end
    ext = '.dat';
    unit_nr = numel(unitnames);
    unit_cell_ids = cell(unit_nr, 1);
    for unit = 1:unit_nr
        centroids = dlmread(strcat(dirc, unitnames{unit}, ext));
        centroids(:,3) = centroids(:,3)*-1;
        disty = pdist2(G.cells.centroids(1:G.layerSize,2), centroids(:,2));
        distz = pdist2(G.cells.centroids(1:G.layerSize,3), centroids(:,3));
        dist = sqrt(disty.^2 + distz.^2);
        
        [~, unit_cell_ids{unit}] = min(dist,[],1);
        unit_cell_ids{unit} = unique(unit_cell_ids{unit});
    end
    
    % reservoir LM2 cells left to Main Fault
    if strcmp(mesh.type, 'sc2')
        if mesh.reduce == 0
            idr = 59:61; ida = 62; idm = 63;
            unit_cell_ids{1} = unique([unit_cell_ids{[1,idr]}]);     % res_LM2 has all cells
            unit_cell_ids{2} = unique([unit_cell_ids{[2,ida]}]);     % amphb    "
            unit_cell_ids{3} = unique([unit_cell_ids{[3,idm]}]);     % MMUM     "
        elseif mesh.ampLyr == 21
            idr = 34:36; ida = 59; idm = 60;
            unit_cell_ids{1} = unique([unit_cell_ids{[1,idr]}]);     % res_LM2 has all cells
            for rightLyrId=1:mesh.ampLyr
                unit_cell_ids{rightLyrId+1} = unique([unit_cell_ids{[rightLyrId+1,rightLyrId+36]}]);
                unitnames{rightLyrId+1} = strcat('A', num2str(rightLyrId));
            end
            unit_cell_ids{58} = unique([unit_cell_ids{[58,ida]}]);    % amphb    "
            unit_cell_ids{23} = unique([unit_cell_ids{[23,idm]}]);   % MMUM     "
        end
            
    elseif strcmp(mesh.type, 'sc21')
        if mesh.reduce == 0
            idr = 54:56; ida = 57; idm = 58;
        elseif mesh.reduce == 1
            idr = 51:53; ida = 54; idm = 55;
        end
        unit_cell_ids{1} = unique([unit_cell_ids{[1,idr]}]);     % res_LM2 has all cells
        unit_cell_ids{2} = unique([unit_cell_ids{[2,ida]}]);     % amphb    "
        unit_cell_ids{3} = unique([unit_cell_ids{[3,idm]}]);     % MMUM     "
    else
        centroidsL(:,3) = centroidsL(:,3)*-1;
        disty = pdist2(G.cells.centroids(1:G.layerSize,2), centroidsL(:,2));
        distz = pdist2(G.cells.centroids(1:G.layerSize,3), centroidsL(:,3));
        dist = sqrt(disty.^2 + distz.^2);
        [~, resLM2Left] = min(dist,[],1);
        resLM2Left = unique(resLM2Left);
    
        % grid refinement area
        centroidsRef(:,3) = centroidsRef(:,3)*-1;
        disty = pdist2(G.cells.centroids(1:G.layerSize,2), centroidsRef(:,2));
        distz = pdist2(G.cells.centroids(1:G.layerSize,3), centroidsRef(:,3));
        dist = sqrt(disty.^2 + distz.^2);
        [~, resLM2ref] = min(dist,[],1);
        resLM2ref = unique(resLM2ref);
        
        % Add cells from left and refined area in case they are not already in
        % resLM2
        if strcmp(mesh.type, 'ref300') || strcmp(mesh.type, 'ref200') || ...
                strcmp(mesh.type, 'ref100') || ...
                all([strcmp(mesh.type, 'ref50'), mesh.reduce == 0]) || ...
                all([strcmp(mesh.type, 'ref40'), mesh.reduce == 0]) || ...
                all([strcmp(mesh.type, 'ref25'), mesh.reduce == 0]) || ...
                all([strcmp(mesh.type, 'ref20'), mesh.reduce == 0]) || ...
                all([strcmp(mesh.type, 'ref15'), mesh.reduce == 0]) || ...
                all([strcmp(mesh.type, 'ref12.5'), mesh.reduce == 0])
            unit_cell_ids{4} = [unit_cell_ids{4} resLM2Left resLM2ref];
        elseif strcmp(mesh.type, 'ref25') || strcmp(mesh.type, 'ref20') || ...
                strcmp(mesh.type, 'ref15') || strcmp(mesh.type, 'ref12.5')
            unit_cell_ids{1} = [unit_cell_ids{1} resLM2Left resLM2ref];
        end
    end
    
    
    
%% toy base porosity and permeability. See permTensor.m
        rock.poro = repmat(.99, [G.cells.num, 1]);
        rock.perm = repmat(300*milli*darcy(), [G.cells.num, 6]);
        rock.perm(:,[2 3 5]) = 0;
    
%     find missing cells by plotting
%     1. Plot base poro and change just one unit at a time, visualize and
%     count number of missed cells.
%     2. Plot cell numbers by area, and assign missed cells
    
%       modify poro, perm differently for each layer
        rock.poro(unit_cell_ids{1}, 1) = .75;
        rock.poro([unit_cell_ids{2:2:22}], 1) = .25;
        rock.poro([unit_cell_ids{3:2:21}], 1) = .65;
        rock.poro(unit_cell_ids{23}, 1) = .1;
        rock.poro(unit_cell_ids{24}, 1) = .8;
        
        rock.poro(unit_cell_ids{25}, 1) = .05;
        rock.poro(unit_cell_ids{26}, 1) = .85;
        rock.poro(unit_cell_ids{27}, 1) = .15;
        rock.poro(unit_cell_ids{28}, 1) = .75;
        rock.poro(unit_cell_ids{29}, 1) = .4;
        rock.poro(unit_cell_ids{30}, 1) = .8;
        rock.poro(unit_cell_ids{31}, 1) = .35;
        rock.poro(unit_cell_ids{32}, 1) = .5;
        rock.poro(unit_cell_ids{33}, 1) = .2;
        
        rock.poro(unit_cell_ids{58}, 1) = .45;
%     
        figure(2); clf, colormap jet
        %plotToolbar(G, rock.poro)
        plotCellData(G, rock.poro, 'EdgeColor', [0.8 0.8 0.8]);
        set(gca,'Xdir','reverse'); axis equal off; view([90 0])
        
        txtargs = {'color', 'k', 'FontSize', 9, 'HorizontalAlignment', 'center'};
        allc = 1:G.layerSize;
        currentc = [unit_cell_ids{1:end}];
        missingc = setdiff(allc,currentc);
        xc = zeros(numel(missingc),1);
        hold on; 
        text(xc, G.cells.centroids(missingc,2), G.cells.centroids(missingc,3), num2str(missingc'), txtargs{:}); hold off
% %     
%         txtargs = {'color', 'k', 'FontSize', 9, 'HorizontalAlignment', 'center'};
%         xv = [0 max(G.cells.centroids(G.layerSize, 1)) + 5]; yv = [20000 35000]; zv = [000 400];
%         idxx = all([G.cells.centroids(:,1)>xv(1), G.cells.centroids(:,1)<xv(2)], 2);
%         idy = all([G.cells.centroids(:,2)>yv(1), G.cells.centroids(:,2)<yv(2)], 2);
%         idz = all([G.cells.centroids(:,3)>zv(1), G.cells.centroids(:,3)<zv(2)], 2);
%         id = all([idxx, idy, idz], 2);
%         cells = 1:G.cells.num; xc = zeros(G.cells.num,1);
%         hold on; 
%         text(xc(id), G.cells.centroids(id,2), G.cells.centroids(id,3), num2str(cells(id)'), txtargs{:}); hold off
%         ylim([yv(1) yv(2)]); zlim([zv(1) zv(2)])
%%

    % assign missing cells in seal layers ALWAYS CHECK FOR REPEATED CELLS
    switch mesh.type
        case 'sc2'
            if mesh.reduce == 1                                
                mis_resLM2 = [4096 4884 12495];
                mis_amp = [19047];
                mis_MMUM = [5592 7688 8868 9237 9250 10453 21928 22092 22285 22326];
                mis_you = [5868 6527 6716];
                mis_al2 = 12806; mis_al3 = 14309; mis_al4 = 14316; mis_al5 = 14325;
                mis_al6 = 14374; mis_al7 = 14410; mis_al8 = 167; mis_al9 = 14438;
                mis_al10 = 14637; mis_al12 = 14671; mis_al13 = 14719; mis_al14 = 14728;
                mis_al16 = 11075;
                mis_ar1  = 13325;

                unit_cell_ids{1} = unique([unit_cell_ids{1} mis_resLM2]);
                unit_cell_ids{58} = unique([unit_cell_ids{58} mis_amp]);
                unit_cell_ids{23} = unique([unit_cell_ids{23} mis_MMUM]);
                unit_cell_ids{24} = unique([unit_cell_ids{24} mis_you]);
                
                unit_cell_ids{2} = unique([unit_cell_ids{2} mis_ar1]);               
                unit_cell_ids{3} = unique([unit_cell_ids{3} mis_al2]); unit_cell_ids{38} = unique([unit_cell_ids{38} mis_al2]);
                unit_cell_ids{4} = unique([unit_cell_ids{4} mis_al3]); unit_cell_ids{39} = unique([unit_cell_ids{39} mis_al3]);
                unit_cell_ids{5} = unique([unit_cell_ids{5} mis_al4]); unit_cell_ids{40} = unique([unit_cell_ids{40} mis_al4]);
                unit_cell_ids{6} = unique([unit_cell_ids{6} mis_al5]); unit_cell_ids{41} = unique([unit_cell_ids{41} mis_al5]);
                unit_cell_ids{7} = unique([unit_cell_ids{7} mis_al6]); unit_cell_ids{42} = unique([unit_cell_ids{42} mis_al6]);
                unit_cell_ids{8} = unique([unit_cell_ids{8} mis_al7]); unit_cell_ids{43} = unique([unit_cell_ids{43} mis_al7]);
                unit_cell_ids{9} = unique([unit_cell_ids{9} mis_al8]); unit_cell_ids{44} = unique([unit_cell_ids{44} mis_al8]);
                unit_cell_ids{10} = unique([unit_cell_ids{10} mis_al9]); unit_cell_ids{45} = unique([unit_cell_ids{45} mis_al9]);
                unit_cell_ids{11} = unique([unit_cell_ids{11} mis_al10]); unit_cell_ids{46} = unique([unit_cell_ids{46} mis_al10]);
                unit_cell_ids{13} = unique([unit_cell_ids{13} mis_al12]); unit_cell_ids{48} = unique([unit_cell_ids{48} mis_al12]);
                unit_cell_ids{14} = unique([unit_cell_ids{14} mis_al13]); unit_cell_ids{49} = unique([unit_cell_ids{49} mis_al13]);
                unit_cell_ids{15} = unique([unit_cell_ids{15} mis_al14]); unit_cell_ids{50} = unique([unit_cell_ids{50} mis_al14]);
                unit_cell_ids{17} = unique([unit_cell_ids{17} mis_al16]); unit_cell_ids{52} = unique([unit_cell_ids{52} mis_al16]);
                
                unit_cell_ids{34} = unique([unit_cell_ids{34} 4096 12495]);   % rl
                unit_cell_ids{35} = unique([unit_cell_ids{35} ]);             % refl
                unit_cell_ids{36} = unique([unit_cell_ids{36} ]);             % ref
                %unit_cell_ids{59} = unique([unit_cell_ids{59} ]);            % amphbleft
                unit_cell_ids{60} = unique([unit_cell_ids{60} 10453]);        % MMUM_left
                
            elseif mesh.reduce == 0
                mis_older  = [19519];
                mis_resLM1 = [842 4456 12905 13013 13073 13282 13335 19118 19658];
                mis_marga  = [1150 18097 18130 ];
                
                mis_resLM2 = [4380 5542 13516];
                mis_amp = [21484];
                mis_MMUM = [6250 8346 9526 9895 9908 11111 24365 24529 24722 24763];
                mis_you = [6526 7185 7374];
                mis_al2 = 13827; mis_al3 = 15367; mis_al4 = 15374; mis_al5 = 15383;
                mis_al6 = 15432; mis_al7 = 15468; mis_al8 = 177; mis_al9 = 15496;
                mis_al10 = 15695; mis_al12 = 15729; mis_al13 = 15777; mis_al14 = 15786;
                mis_al16 = 11733;
                mis_ar1  = 14367;
                %unit_cell_ids{10}(unit_cell_ids{10}==12912)=[];
                
                unit_cell_ids{56} = unique([unit_cell_ids{56} mis_older]);
                unit_cell_ids{57} = unique([unit_cell_ids{57} mis_resLM1]);
                unit_cell_ids{58} = unique([unit_cell_ids{58} mis_marga]);
                unit_cell_ids{1} = unique([unit_cell_ids{1} mis_resLM2]);
                unit_cell_ids{2} = unique([unit_cell_ids{2} mis_amp]);
                unit_cell_ids{3} = unique([unit_cell_ids{3} mis_MMUM]);
                unit_cell_ids{4} = unique([unit_cell_ids{4} mis_you]);
                
                unit_cell_ids{14} = unique([unit_cell_ids{14} mis_ar1]);
                
                unit_cell_ids{36} = unique([unit_cell_ids{36} mis_al2]);
                unit_cell_ids{37} = unique([unit_cell_ids{37} mis_al3]);
                unit_cell_ids{38} = unique([unit_cell_ids{38} mis_al4]);
                unit_cell_ids{39} = unique([unit_cell_ids{39} mis_al5]);
                unit_cell_ids{40} = unique([unit_cell_ids{40} mis_al6]);
                unit_cell_ids{41} = unique([unit_cell_ids{41} mis_al7]);
                unit_cell_ids{42} = unique([unit_cell_ids{42} mis_al8]);
                unit_cell_ids{43} = unique([unit_cell_ids{43} mis_al9]);
                unit_cell_ids{44} = unique([unit_cell_ids{44} mis_al10]);
                unit_cell_ids{46} = unique([unit_cell_ids{46} mis_al12]);
                unit_cell_ids{47} = unique([unit_cell_ids{47} mis_al13]);
                unit_cell_ids{48} = unique([unit_cell_ids{48} mis_al14]);
                unit_cell_ids{50} = unique([unit_cell_ids{50} mis_al16]);
                
                unit_cell_ids{59} = unique([unit_cell_ids{59} 4380 13516]);   % rl
                unit_cell_ids{61} = unique([unit_cell_ids{61} ]);             % refl
                unit_cell_ids{60} = unique([unit_cell_ids{60} ]);             % ref
                %unit_cell_ids{57} = unique([unit_cell_ids{62} ]);            % amphbleft
                unit_cell_ids{63} = unique([unit_cell_ids{63} 11111]);        % MMUM_left
            end
            
        case 'sc21'
            if mesh.reduce == 1                
                mis_resLM2 = [919 1308 14077 14226 14417 16377 17822 21366];
                mis_amp = [6186];
                mis_MMUM = [2192 5677 9430 10687 15925];
                mis_you = [6796 7086];
                mis_al2 = 3430; mis_al3 = 13935; mis_al4 = 13942; mis_al5 = 3824;
                mis_al6 = 13741; mis_al7 = 15296; mis_al8 = 15370; mis_al9 = 3855;
                mis_al10 = 3858; mis_al12 = 15457; mis_al13 = 15495; mis_al14 = 15692;
                mis_al16 = 3941;
                
                unit_cell_ids{1} = unique([unit_cell_ids{1} mis_resLM2]);
                unit_cell_ids{2} = unique([unit_cell_ids{2} mis_amp]);
                unit_cell_ids{3} = unique([unit_cell_ids{3} mis_MMUM]);
                unit_cell_ids{4} = unique([unit_cell_ids{4} mis_you]);
                
                unit_cell_ids{31} = unique([unit_cell_ids{31} mis_al2]);
                unit_cell_ids{32} = unique([unit_cell_ids{32} mis_al3]);
                unit_cell_ids{33} = unique([unit_cell_ids{33} mis_al4]);
                unit_cell_ids{34} = unique([unit_cell_ids{34} mis_al5]);
                unit_cell_ids{35} = unique([unit_cell_ids{35} mis_al6]);
                unit_cell_ids{36} = unique([unit_cell_ids{36} mis_al7]);
                unit_cell_ids{37} = unique([unit_cell_ids{37} mis_al8]);
                unit_cell_ids{38} = unique([unit_cell_ids{38} mis_al9]);
                unit_cell_ids{39} = unique([unit_cell_ids{39} mis_al10]);
                unit_cell_ids{41} = unique([unit_cell_ids{41} mis_al12]);
                unit_cell_ids{42} = unique([unit_cell_ids{42} mis_al13]);
                unit_cell_ids{43} = unique([unit_cell_ids{43} mis_al14]);
                unit_cell_ids{45} = unique([unit_cell_ids{45} mis_al16]);
                
                unit_cell_ids{51} = unique([unit_cell_ids{51} 14417 16377 17822]);       % rl
                unit_cell_ids{53} = unique([unit_cell_ids{53} 14077 14226]);             % refl
                unit_cell_ids{52} = unique([unit_cell_ids{52} 919 21366]);               % ref
                %unit_cell_ids{54} = unique([unit_cell_ids{54} ]);                        % amphbleft
                unit_cell_ids{55} = unique([unit_cell_ids{55} 9430 ]);
                
            elseif mesh.reduce == 0
                mis_older  = [5171];
                mis_resLM1 = [13885 13888 14129 14134 14223 14225 3545 19141 5085 20597];
                mis_marga  = [300 19168 19304];
                
                mis_resLM2 = [988 1480 15103 15265 15456 17421 18866 23803 ];
                mis_amp = [6822];
                mis_MMUM = [2364 6313 10066 11323 16969];
                mis_you = [7432 7722];
                mis_al2 = 3673; mis_al3 = 14961; mis_al4 = 14968; mis_al5 = 4080;
                mis_al6 = 14767; mis_al7 = 16340; mis_al8 = 16414; mis_al9 = 4111;
                mis_al10 = 4114; mis_al12 = 16501; mis_al13 = 16539; mis_al14 = 16736;
                mis_al16 = 4197;
                %unit_cell_ids{10}(unit_cell_ids{10}==12912)=[];
                
                unit_cell_ids{51} = unique([unit_cell_ids{51} mis_older]);
                unit_cell_ids{52} = unique([unit_cell_ids{52} mis_resLM1]);
                unit_cell_ids{53} = unique([unit_cell_ids{53} mis_marga]);
                unit_cell_ids{1} = unique([unit_cell_ids{1} mis_resLM2]);
                unit_cell_ids{2} = unique([unit_cell_ids{2} mis_amp]);
                unit_cell_ids{3} = unique([unit_cell_ids{3} mis_MMUM]);
                unit_cell_ids{4} = unique([unit_cell_ids{4} mis_you]);
                
                unit_cell_ids{31} = unique([unit_cell_ids{31} mis_al2]);
                unit_cell_ids{32} = unique([unit_cell_ids{32} mis_al3]);
                unit_cell_ids{33} = unique([unit_cell_ids{33} mis_al4]);
                unit_cell_ids{34} = unique([unit_cell_ids{34} mis_al5]);
                unit_cell_ids{35} = unique([unit_cell_ids{35} mis_al6]);
                unit_cell_ids{36} = unique([unit_cell_ids{36} mis_al7]);
                unit_cell_ids{37} = unique([unit_cell_ids{37} mis_al8]);
                unit_cell_ids{38} = unique([unit_cell_ids{38} mis_al9]);
                unit_cell_ids{39} = unique([unit_cell_ids{39} mis_al10]);
                unit_cell_ids{41} = unique([unit_cell_ids{41} mis_al12]);
                unit_cell_ids{42} = unique([unit_cell_ids{42} mis_al13]);
                unit_cell_ids{43} = unique([unit_cell_ids{43} mis_al14]);
                unit_cell_ids{45} = unique([unit_cell_ids{45} mis_al16]);
                
                unit_cell_ids{54} = unique([unit_cell_ids{54} 15456 17421 18866]);       % rl
                unit_cell_ids{56} = unique([unit_cell_ids{56} 15103 15265]);             % refl
                unit_cell_ids{55} = unique([unit_cell_ids{55} 988 23803]);               % ref
                %unit_cell_ids{57} = unique([unit_cell_ids{58} ]);                        % amphbleft
                unit_cell_ids{58} = unique([unit_cell_ids{58} 10066]);                        % MMUM_left
            end    
            
        case 'ref12.5_upper'
            mis_resLM2 = [];
            mis_amp = [];   
            mis_fresLM2  = [];
            mis_fresLM2_amp  = [];
            unit_cell_ids{5}(unit_cell_ids{5}==[])=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_resLM2]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_amp]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_fresLM2]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left ]);
            resLM2ref  = [resLM2ref];
        case 'ref12.5'
            mis_older = [1421 21401];
            mis_resLM1 = [5323 22732];
            mis_marg = [21454 22059];
            mis_resLM2 = [68 3725 7387 15014 15035 19980 21121 21355];
            mis_amp = [2417 10119 10145 12777 28640 28989 29481];        
            mis_MMUM = [8568 9378 9714 10228];
            mis_you = [508 30458];
            mis_fresLM2  = [17161];
            mis_fresLM2_amp  = [12910 12912];
            unit_cell_ids{10}(unit_cell_ids{10}==12912)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_older]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_resLM1]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_marg]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_resLM2]);
            unit_cell_ids{5} = unique([unit_cell_ids{5} mis_amp]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_MMUM]);
            unit_cell_ids{7} = unique([unit_cell_ids{7} mis_you]);
            unit_cell_ids{8} = unique([unit_cell_ids{8} mis_fresLM2]);
            unit_cell_ids{9} = unique([unit_cell_ids{9} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left 3725 15014 15035]);
            resLM2ref  = [resLM2ref];
        case 'ref15_upper'
            mis_resLM2 = [7272 9771];
            mis_amp = [1244 7004 18322];        
            %mis_MMUM = [];
            %mis_you = [];
            mis_fresLM2_amp  = [1723 7267];
            unit_cell_ids{5}(unit_cell_ids{5}==1723)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_resLM2]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_amp]);
            %unit_cell_ids{3} = unique([unit_cell_ids{3} mis_MMUM]);
            %unit_cell_ids{4} = unique([unit_cell_ids{4} mis_you]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left 7272 9771]);
            resLM2ref  = [resLM2ref];
        case 'ref15'
            mis_older = [18756];
            mis_resLM1 = [1151 17247];
            mis_marg = [17402 17577];
            mis_resLM2 = [10836 13693];
            mis_amp = [2032 10568 23563];        
            mis_MMUM = [1948 7476 7775 8476];
            mis_you = [6762 24858];
            mis_fresLM2_amp  = [2574 10831];
            unit_cell_ids{10}(unit_cell_ids{10}==2574)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_older]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_resLM1]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_marg]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_resLM2]);
            unit_cell_ids{5} = unique([unit_cell_ids{5} mis_amp]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_MMUM]);
            unit_cell_ids{7} = unique([unit_cell_ids{7} mis_you]);
            unit_cell_ids{9} = unique([unit_cell_ids{9} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left 10836 13693]);
            resLM2ref  = [resLM2ref];
        case 'ref20_upper'
            mis_resLM2 = [];
            mis_amp = [6182 7257 7751];        
            mis_MMUM = [5442 5470 7348 7502];
            mis_you = [4706 16679];
            mis_fresLM2_amp  = [6086 6076];
            unit_cell_ids{7}(unit_cell_ids{7}==6076)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_resLM2]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_amp]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_MMUM]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_you]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left ]);
            resLM2ref  = [resLM2ref];
        case 'ref20'
            mis_older = [3573];
            mis_resLM1 = [207 862];
            mis_marg = [9906 10843 13242 13732];
            mis_resLM2 = [];
            mis_amp = [6756 7831 8325];        
            mis_MMUM = [6016 6044 7922 8076];
            mis_you = [5280 18889];
            mis_fresLM2_amp  = [6660 6650];
            unit_cell_ids{10}(unit_cell_ids{10}==6650)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_older]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_resLM1]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_marg]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_resLM2]);
            unit_cell_ids{5} = unique([unit_cell_ids{5} mis_amp]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_MMUM]);
            unit_cell_ids{7} = unique([unit_cell_ids{7} mis_you]);
            unit_cell_ids{9} = unique([unit_cell_ids{9} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left ]);
            resLM2ref  = [resLM2ref];
        case 'ref25_upper'
            mis_resLM2 = [8310 12198];
            mis_amp = [5881 5227];        
            mis_MMUM = [247 4704 6188];
            mis_you = [3923 14301];
            mis_fresLM2_amp  = [5797 5801];
            unit_cell_ids{7}(unit_cell_ids{7}==5797)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_resLM2]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_amp]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_MMUM]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_you]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left 8310]);
            resLM2ref  = [resLM2ref];
        case 'ref25'
            mis_older = [12323];
            mis_resLM1 = [11379 12294];
            mis_marg = [11293 11784];
            mis_resLM2 = [9282 14403];
            mis_amp = [6490 5836];        
            mis_MMUM = [278 5313 6797];
            mis_you = [4532 16506];
            mis_fresLM2_amp  = [6410 6406];
            unit_cell_ids{10}(unit_cell_ids{10}==6406)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_older]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_resLM1]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_marg]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_resLM2]);
            unit_cell_ids{5} = unique([unit_cell_ids{5} mis_amp]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_MMUM]);
            unit_cell_ids{7} = unique([unit_cell_ids{7} mis_you]);
            unit_cell_ids{9} = unique([unit_cell_ids{9} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left 9282]);
            resLM2ref  = [resLM2ref];
        case 'ref40'
            mis_older = [2487];
            mis_resLM1 = [9957 7488];
            mis_marg = [9171 9101];
            mis_resLM2 = [1726 6956];
            mis_amp = [1233 1547 5848 8170 12521];        
            mis_MMUM = [1089 4335 4853 5237 12555];
            mis_you = [3656 13400];
            mis_fresLM2_amp  = [6027 6146];
            unit_cell_ids{10}(unit_cell_ids{10}==6027)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_older]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_resLM1]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_marg]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_resLM2]);
            unit_cell_ids{5} = unique([unit_cell_ids{5} mis_amp]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_MMUM]);
            unit_cell_ids{7} = unique([unit_cell_ids{7} mis_you]);
            unit_cell_ids{9} = unique([unit_cell_ids{9} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left 1726 6956]);
            resLM2ref  = [resLM2ref ];
        case 'ref50'
            mis_older = [9510];
            mis_resLM1 = [7001 9457];
            mis_marg = [8462 8595];
            mis_resLM2 = [1707 6729 7646 11357];
            mis_amp = [1922 4829 5682 11614];        
            mis_MMUM = [922 4131 4158 4386];
            mis_you = [3506 12729];
            mis_fresLM2_amp  = [1429 6054];
            unit_cell_ids{10}(unit_cell_ids{10}==1429)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_older]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_resLM1]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_marg]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_resLM2]);
            unit_cell_ids{5} = unique([unit_cell_ids{5} mis_amp]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_MMUM]);
            unit_cell_ids{7} = unique([unit_cell_ids{7} mis_you]);
            unit_cell_ids{9} = unique([unit_cell_ids{9} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left 1707 6729]);
            resLM2ref  = [resLM2ref 7646];
        case 'ref100'
            mis_older = [8860];
            mis_resLM1 = [8675 414];
            mis_marg = [6683 8031];
            mis_resLM2 = [1758 6312 7083 7181];
            mis_amp = [5241 5308 10958];        
            mis_MMUM = [3848 4123 4161];
            mis_you = [3247 11422];
            mis_fresLM2_amp  = [1772 5519];
            unit_cell_ids{10}(unit_cell_ids{10}==5519)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_older]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_resLM1]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_marg]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_resLM2]);
            unit_cell_ids{5} = unique([unit_cell_ids{5} mis_amp]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_MMUM]);
            unit_cell_ids{7} = unique([unit_cell_ids{7} mis_you]);
            unit_cell_ids{9} = unique([unit_cell_ids{9} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left 6312 7083]);
            resLM2ref  = [resLM2ref 1758 7181];
        case 'ref200'
            mis_older = [6833];
            mis_resLM1 = [6233 6555];
            mis_marg = [6032];
            mis_resLM2 = [5659 7098 7952];
            mis_amp = [73 4779];        
            mis_MMUM = [3213 3597 4071];
            mis_you = [2148 2456];
            mis_fresLM2_amp  = [5142 5183];
            unit_cell_ids{10}(unit_cell_ids{10}==5183)=[];
            
            unit_cell_ids{1} = unique([unit_cell_ids{1} mis_older]);
            unit_cell_ids{2} = unique([unit_cell_ids{2} mis_resLM1]);
            unit_cell_ids{3} = unique([unit_cell_ids{3} mis_marg]);
            unit_cell_ids{4} = unique([unit_cell_ids{4} mis_resLM2]);
            unit_cell_ids{5} = unique([unit_cell_ids{5} mis_amp]);
            unit_cell_ids{6} = unique([unit_cell_ids{6} mis_MMUM]);
            unit_cell_ids{7} = unique([unit_cell_ids{7} mis_you]);
            unit_cell_ids{9} = unique([unit_cell_ids{9} mis_fresLM2_amp]);
            
            resLM2Left = unique([resLM2Left 5659 7952]);
            %resLM2ref  = [resLM2ref ];
        case 'ref300'
            mis_older = [5707 6240];
            mis_resLM1 = [5858 1471];
            mis_marg = [341 331];
            mis_resLM2 = [1802 5190 7170 7199];
            mis_amp = [1121 4397 7362];        
            mis_MMUM = [2941 4773 5079];
            mis_you = [2269 7883];
            mis_fresLM2_amp  = [7341 7342];
            unit_cell_ids{10}(unit_cell_ids{10}==7342)=[];
            
            unit_cell_ids{1} = [unit_cell_ids{1} mis_older];
            unit_cell_ids{2} = [unit_cell_ids{2} mis_resLM1];
            unit_cell_ids{3} = [unit_cell_ids{3} mis_marg];
            unit_cell_ids{4} = [unit_cell_ids{4} mis_resLM2];
            unit_cell_ids{5} = [unit_cell_ids{5} mis_amp];
            unit_cell_ids{6} = [unit_cell_ids{6} mis_MMUM];
            unit_cell_ids{7} = [unit_cell_ids{7} mis_you];
            unit_cell_ids{9} = unique([unit_cell_ids{9} mis_fresLM2_amp]);
            
            resLM2Left = [resLM2Left 5190 7170];
            resLM2ref  = [resLM2ref 1802 7199];
            
        case 'coarse'
            mis_older = [8889];
            mis_resLM1 = [6477 8668];
            mis_marg = [409 8088];
            mis_resLM2 = [436 7188 7084];
            mis_amp = [2523 5643 7291 10908 10856 11076];        
            mis_MMUM = [1134 3870 4192];
            mis_you = [816 11874];
            mis_fresLM2_amp  = 5604;
            unit_cell_ids{10}(unit_cell_ids{10}==5598)=[];
            
            unit_cell_ids{1} = [unit_cell_ids{1} mis_older];
            unit_cell_ids{2} = [unit_cell_ids{2} mis_resLM1];
            unit_cell_ids{3} = [unit_cell_ids{3} mis_marg];
            unit_cell_ids{4} = [unit_cell_ids{4} mis_resLM2];
            unit_cell_ids{5} = [unit_cell_ids{5} mis_amp];
            unit_cell_ids{6} = [unit_cell_ids{6} mis_MMUM];
            unit_cell_ids{7} = [unit_cell_ids{7} mis_you];
            unit_cell_ids{9} = [unit_cell_ids{9} mis_fresLM2_amp];
            
            resLM2Left = [resLM2Left 436 7084];
    
        case 'coarseUpper'
            mis_resLM2 = [362 6123 6019];
            mis_amp = [1950 5070 6226 8707 8655 8875];        
            mis_MMUM = [970 3297 3619];
            mis_you = [652 9673];
            mis_fresLM2_amp  = 5031;
            unit_cell_ids{7}(unit_cell_ids{7}==5025)=[];
            
            unit_cell_ids{1} = [unit_cell_ids{1} mis_resLM2];
            unit_cell_ids{2} = [unit_cell_ids{2} mis_amp];
            unit_cell_ids{3} = [unit_cell_ids{3} mis_MMUM];
            unit_cell_ids{4} = [unit_cell_ids{4} mis_you];
            unit_cell_ids{6} = [unit_cell_ids{6} mis_fresLM2_amp];
            
            resLM2Left = [resLM2Left 362 6019];
        case 'fine2'
            % later
    end

    % complete cell assignment to extruded dimension
    id_add = 1:G.numLayers-1;
    ucid = unit_cell_ids;
    for unit = 1:unit_nr
        for lyr = 1:numel(id_add)
        unit_cell_ids{unit} = [unit_cell_ids{unit} ucid{unit}+G.layerSize*id_add(lyr)];
        end
    end
    
%     lyr1 = resLM2Left;
%     lyr2 = resLM2ref;
%     lyr4 = elemOnTop;  
%     for lyr = 1:numel(id_add)
%         resLM2Left = [resLM2Left lyr1+G.layerSize*id_add(lyr)];
%         resLM2ref = [resLM2ref lyr2+G.layerSize*id_add(lyr)];
%         elemOnTop = [elemOnTop lyr6+G.layerSize*id_add(lyr)];
%     end
    
    resLM2Left = unit_cell_ids{34}; 
    resLM2ref = unit_cell_ids{35}; 
    resLM2refl = unit_cell_ids{36}; 
    amphLeft = unit_cell_ids{59}; 
    MMUMleft = unit_cell_ids{60}; 
    
    
    % add full fault for convenience (plotting, etc)
    unitnames{end + 1} = 'fault_all';
    if mesh.reduce == 2
        unit_cell_ids{end+1} = [unit_cell_ids{end-2} unit_cell_ids{end-1} unit_cell_ids{end}];
    else
        unit_cell_ids{end+1} = [unit_cell_ids{end-5} unit_cell_ids{end-4} unit_cell_ids{end-3} ...
                                unit_cell_ids{end-2} unit_cell_ids{end-1} unit_cell_ids{end}];
    end
    if strcmp(mesh.type, 'sc21')
        unit_cell_ids{end+1} = [unit_cell_ids{5:8}];
    elseif strcmp(mesh.type, 'sc2')
        unit_cell_ids{end+1} = [unit_cell_ids{25:33}];
    end

%     % check individual unit
%     clf
%     figure(2)
%     plotCellData(G, rock.perm(:,1), [unit_cell_ids{8} unit_cell_ids{9} ...
%                  unit_cell_ids{10} unit_cell_ids{11} unit_cell_ids{12} ...
%                  unit_cell_ids{13}], 'EdgeColor', [0.8 0.8 0.8]);
%     plotCellData(G, rock.poro, unit_cell_ids{1}, 'EdgeColor', [0.8 0.8 0.8]);
%     axis equal
%     set(gca,'Xdir','reverse'); axis off; view([90 0])
    
    % save data
    disp('WARNING: You are about to save unit_cell_ids and resLM2Left data. Press any key to continue')
    pause
    fid = ['ucids_' mesh.type '.mat'];
    save(fid,'unit_cell_ids','unitnames');
    save('ucids_resLM2Left.mat', 'resLM2Left')
    save('ucids_resLM2ref.mat', 'resLM2ref')
    save('ucids_resLM2refl.mat', 'resLM2refl')
    save('ucids_amphLeft.mat', 'amphLeft')
    save('ucids_MMUMleft.mat', 'MMUMleft')
    return
end

if obtained.fnodeset == 0
   if strcmp(mesh.type,'coarse') || strcmp(mesh.type,'coarseUpper')
       %meshinfo =
       %ncinfo('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/extr_tri/gom_mz_tri_fdiv.exo'); % to know which nodeset
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gom_mz_tri_fdiv.exo','node_ns5');
   elseif strcmp(mesh.type, 'ref300')
       %ncinfo('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref300.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref300.exo','node_ns5');
   elseif strcmp(mesh.type, 'ref200')
       %ncinfo('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref200.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref200.exo','node_ns5');
   elseif strcmp(mesh.type, 'ref100')
       %ncinfo('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref100.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref100.exo','node_ns5');
   elseif strcmp(mesh.type, 'ref50')
       %ncinfo('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref50.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref50.exo','node_ns5');
   elseif strcmp(mesh.type, 'ref40')
       %ncinfo('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref40.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref40.exo','node_ns5');
   elseif strcmp(mesh.type, 'ref25')
       %ncinfo('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref25.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref25.exo','node_ns5');
   elseif strcmp(mesh.type, 'ref20')
       %ncinfo('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref25.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref20.exo','node_ns5');
   elseif strcmp(mesh.type, 'ref15')
       %ncinfo('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref25.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref15.exo','node_ns5');
   elseif strcmp(mesh.type, 'ref12.5')
       %ncinfo('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref25.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/revised1_ref/gom_mz_tri_fdiv_ref12.5.exo','node_ns5');
   elseif strcmp(mesh.type, 'sc21')
       %ncinfo('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/scenario2/gom_mz_tri_s21.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/scenario2/gom_mz_tri_s21.exo','node_ns5');
   elseif strcmp(mesh.type, 'sc21')
       %ncinfo('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/scenario2/gom_mz_tri_s2_throwDiv.exo');
       fnod = ncread('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/meshes/scenario2/gom_mz_tri_s2_throwDiv.exo','node_ns5');
   else
       error('get nodes for fine mesh!')
   end
   fnodcoord = G.nodes.coords(fnod, :);
   % save data
   disp('WARNING: You are about to save fnodcoord data. Press any key to continue')
   pause
   fid = 'fnodcoord.mat';
   save(fid, 'fnodcoord');
end

if isfield(mesh, 'ref') && mesh.ref == 1 && obtained.elemOnTop == 0
    fid = 'elemOnTopAreaInX.mat';
    save(fid, 'elemOnTop'); 
end

end
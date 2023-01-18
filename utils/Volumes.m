classdef Volumes
   methods (Static)
       
       function [v_hybrid, v_fine] = discretizationVol(d, model_hybrid, model_fine, state_hybrid_fs, state_fine)
            % Compute CO2 volume inside given discretization region.
            % Args:
            %   d: discretization region
            %   model_hybrid: hybrid model
            %   model_fine: full-dimensional model
            %   state_hybrid_fs: hyhbrid sf reconstructed for fine cells
            %   state_fine: sf from full-dimensional model
            % Returns:
            %   v_hybrid: CO2 volume in discetization region from hybrid model
            %   v_fine: CO2 volume in disretization region from full-dim model
            discr = model_hybrid.G.cells.discretization == d;
            p = model_hybrid.G.partition;
            discr = discr(p);

            snh = state_hybrid_fs.s(:,2);    
            snh = snh(discr);
            snf = state_fine.s(:,2);
            snf = snf(discr);
            pv = poreVolume(model_fine.G, model_fine.rock);
            pv = pv(discr);

            v_hybrid = pv.*snh;
            v_fine = pv.*snf;
       end
            
        
       function [vol_hybrid, vol_fine] = plotSealingLayersVols(model_h, model_f, states_hf, states_f, fig_nr)
            op_h = model_h.operators;
            n = op_h.N;
            discr = model_h.G.cells.discretization;
           
            vol_hybrid = zeros(numel(states_hf), 1);
            vol_fine = zeros(numel(states_f), 1);

            cB = states_hf{1}.cBottom; % bottom cells of semi-perm layers
            veB = states_hf{1}.fBottom; % bottom faces of semi-perm-layers
            veBottom = ismember(n, n(op_h.connections.veToFineVertical | op_h.connections.veTransitionVerticalConn, :), 'rows');
            veB_global = find(veBottom);
            veB_global = veB_global(veB); % indicator for bottom faces for semi-perm layers
            db = unique(discr(cB)); % discretization region for bottom cells of semi-perm layers
            for d=1:numel(db)
                dc = discr(cB) == db(d);
                bottom_conn = veB_global(dc);
                bottom_flux = any(states_hf{end}.vGsum_sMax(bottom_conn));
                if ~bottom_flux
                   continue; % skip discretization region if no CO2 originates from bottom
                end

                fig_nr = fig_nr + 1;
                figure(fig_nr);
                for i=1:numel(states_f)
                    [vol_h, vol_f] = Volumes.discretizationVol(db(d), model_h, model_f, ...
                                                                states_hf{i}, states_f{i});
                    vol_hybrid(i) = sum(vol_h);
                    vol_fine(i) = sum(vol_f);
                end

                plot(1:numel(states_f), vol_fine, 'DisplayName', 'fine')
                hold on
                plot(1:numel(states_hf), vol_hybrid, 'DisplayName', 'hybrid')
                xlabel('Time step');
                ylabel('m^3');
                title(sprintf('CO2 volume in discretization region %d (bottom fluxes)', db(d)));
                legend('location', 'northwest')
            end

            % Plot CO2 volume in discretization regions adjacent to semi-perm layers
            % with horizontal and bottom fluxes
            cH = states_hf{1}.cHorz;
            veH = states_hf{1}.fHorz;
            veHorz = ismember(n, n(op_h.connections.veTransitionHorizontalConn, :), 'rows');
            veH_global = find(veHorz);
            veH_global = veH_global(veH);

            cBH = states_hf{1}.cBottomHorz;
            veBH = states_hf{1}.fBottomHorz;
            veBottomHorz = veBottom | veHorz;
            veBH_global = find(veBottomHorz);
            veBH_global = veBH_global(veBH);

            cBH_all = [cH; cBH];
            veBH_all = [veH_global; veBH_global]; 

            dbh = unique(discr(cBH_all));
            for d=1:numel(dbh)
                dc = discr(cBH_all) == dbh(d);
                bh_conn = veBH_all(dc);
                bh_flux = any(states_hf{end}.vGsum_sMax(bh_conn));
                if ~bh_flux
                   continue; % skip discretization region if no CO2 originates from bottom
                end

                fig_nr = fig_nr + 1;
                figure(fig_nr);
                for i=1:numel(states_f)
                    [vol_h, vol_f] = Volumes.discretizationVol(dbh(d), model_h, model_f, ...
                                                                states_hf{i}, states_f{i});
                    vol_hybrid(i) = sum(vol_h);
                    vol_fine(i) = sum(vol_f);
                end

                plot(1:numel(states_f), vol_fine, 'DisplayName', 'fine')
                hold on
                plot(1:numel(states_hf), vol_hybrid, 'DisplayName', 'hybrid')
                xlabel('Time step');
                ylabel('m^3');
                title(sprintf('CO2 volume in discretization region %d (horizontal fluxes)', dbh(d)));
                legend('location', 'northwest')
            end
       end
       
       function [mean_fh, var_fh] = diffCO2SatCoarse(model_h, model_f, sn_h, sn_f, fig_nr, plot_dir)
       %Compare CO2 saturation between fine and hybrid models 
       %for each coarse region of hybrid model.
       %
       % PARAMETERS:
       %    model_h     - Hybrid model
       %    model_f     - Fine model
       %    sn_h        - co2 saturation in (coarse) cells of hybrid model
       %    sn_f        - co2 saturation in cells of fine model
       %    fig_nr      - figure number for plotting
       %    plot_dir    - directory to save figures
       %
       % RETURNS:
       %    mean_fh     - mean difference in co2 saturation
       %    var_fh      - variance in co2 saturation
            p = model_h.G.partition;

            [vol_unique, ~, idx_sort] = uniquetol(model_h.G.cells.volumes);

            pvh = poreVolume(model_h.G, model_h.rock);
            pvf = poreVolume(model_f.G, model_f.rock);

            sn_f_net = accumarray(p, sn_f.*pvf);
            sn_f_net = sn_f_net ./ pvh;
            
            mean_fh = zeros(numel(vol_unique), 1);
            var_fh = zeros(size(mean_fh));
            % Volume mismatch sorted after increasing volume
            for i=1:numel(vol_unique)    
                iv = idx_sort == i;
                mean_fh(i) = mean(abs(sn_h(iv) - sn_f_net(iv)));    
                var_fh(i) = var(abs(sn_h(iv) - sn_f_net(iv)));
            end

            fig = figure(fig_nr);
            plot(1:numel(mean_fh), mean_fh, 'b', 'DisplayName', 'Mean')
            hold on
            plot(1:numel(var_fh), var_fh, '-r', 'DisplayName', 'Variance')
            % SET VOLUME TICKS !
            xlabel('Coarse cells (increasing volume)')
            title('Absolute difference in CO2 saturation: coarse')
            legend();
            drawnow;
            saveas(fig, strcat(plot_dir, 'diff_sat_coarse'), 'png')
       end
        
       function [diff_fh, var_fh] = diffCO2SatFine(model_h, sn_hf, sn_f, fig_nr, plot_dir)
            discr = model_h.G.cells.discretization;
            discr_u = unique(discr);
            p = model_h.G.partition;

            diff_fh = zeros(numel(discr_u), 1);
            var_fh = zeros(size(diff_fh));
            for i=1:numel(discr_u) % discretization starts at 1 for fine cells
                d = discr_u(i);
                dp = discr(p);
                iv = dp == d;
                diff_fh(i) = mean(abs(sn_hf(iv) - sn_f(iv)));    
                var_fh(i) = var(abs(sn_hf(iv) - sn_f(iv)));
            end

            fig = figure(fig_nr);
            plot(1:numel(diff_fh), diff_fh, 'b', 'DisplayName', 'Mean')
            hold on
            plot(1:numel(var_fh), var_fh, '-r', 'DisplayName', 'Variance')
            % SET VOLUME TICKS !
            xticklabels(discr_u);
            xlabel('Discretization region')
            title('Absolute difference in CO2 saturation: reconstruction')
            legend();
            drawnow;
            saveas(fig, strcat(plot_dir, 'diff_sat_fine'), 'png') 
       end  
       
       function [Vh_exit, Vf_exit] = plotExitedVolumes(model_f, states_hybrid_f, states_fine, ...
                                                        schedule, fig_nr, plot_dir)
            fig = figure(fig_nr);
            dt = schedule.step.val;            
            
            rate = [];
            for i=1:numel(schedule.control)
                num_W = nnz(schedule.step.control == i); % account for all wells
                rate_i = schedule.control(i).W.val * schedule.control(i).W.status; % only add rate if well is on
                rate = cat(1, rate, repmat(rate_i, num_W, 1));
            end
            
            V_inj = cumsum(rate.*dt);
            pvf = poreVolume(model_f.G, model_f.rock);
            
            Vf_net = [];
            Vh_net = [];
            for i=1:numel(states_fine)
               Vf_net = cat(1, Vf_net, sum(states_fine{i}.s(:,2).*pvf)); 
               Vh_net = cat(1, Vh_net, sum(states_hybrid_f{i}.s(:,2).*pvf));
            end

            Vf_exit = V_inj - Vf_net;
            Vh_exit = V_inj - Vh_net;

            plot(1:numel(states_fine), Vf_exit, 'b', 'DisplayName', 'Fine model')
            hold on
            plot(1:numel(states_hybrid_f), Vh_exit, 'r', 'DisplayName', 'Hybrid model')
            xlabel('Time step')
            ylabel('m^3')
            title('CO2 volume exited domain')
            legend('location', 'northwest');
            drawnow;
            saveas(fig, strcat(plot_dir, 'volume_exit'), 'png')
       end
   end
end
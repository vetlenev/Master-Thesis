function [c_VE_Not, c_VE, c_VE_Horz, R] = getResidualFilledCells_test(model, Sn, vG, Sn0, sGnve0)
    % Get VE columns that are filled with CO2 from bottom.
    % Used for correct computation of fluxes in the residual framework.
    %
    % Inputs:
    %   T/B: top/bottom of virtual cells
    %   t/b: top/bottom of column
    %   h: global height of plume in column
    %   h_max: global max depth CO2 has reached
    %   swr/snr: residual water/co2 saturation
    %   veBottomNonzeroSn: boolean array indicating cells in vertical
    %                      ve transition zone, where ve is upper part and ve/fine cells lower part.
    %   Sn: Co2 saturation from COARSE grid    
    
    CG = model.G;
    op = model.operators;    
    isFine = CG.cells.discretization == 1;      
    
    p = CG.partition;
    %Sn_p = Sn(p);
    %snMax_p = snMax(p);
    
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);
    %[a_M, a_R, sG] = getGasSatFromHeight(TT, tt, BB, bb, h, h_max, swr, snr);
    % ------------------- 
      
    % Modify for ve bottom transition
    veTransition = op.connections.veToFineVertical | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;
    veInternalConn = op.connections.veInternalConn;
    veHorzConn = op.connections.veTransitionHorizontalConn; % transitio between two VE cols of different discretizations
    
    c_VE_Not = {};
    c_VE = {};
    c_VE_Horz = {};
    R = {};
    
    n = op.N; 
    %n = CG.faces.neighbors;
    cn = op.N(veTransition, :);
    %c_vihc = op.N(veInternalConn | veHorzConn, :); % candidates for horizontal fluxes
    c_vic = op.N(veInternalConn, :); % make it first work with internal ve before proceeding with ve horizontal conns
    
    % ---
    vG = abs(vG);
    % ---
    
    vG_veTrans = vG(ismember(n, cn, 'rows'));
    % NB: Get ONE horizontal flux value per cn-cell
    %vG_veHorz = vG(ismember(n, c_vihc, 'rows'));   
    vG_veHorz = zeros(numel(vG_veTrans), 1);
    
    for idx=1:2 % CONSIDER REMOVING LOOP: veBottomFilled and veTopEmpty will never both occur as it stand now!
        c = cn(:, idx);    
        isVE_c = ~isFine(c);
           
        if any(isVE_c) % reconstruct specifically for ve col transitioning to a lower layer
            t = CG.cells.topDepth(c);
            T = op.connections.faceTopDepth(veTransition, idx);
            b = CG.cells.bottomDepth(c);
            B = op.connections.faceBottomDepth(veTransition, idx);
            H = b - t; 
            
            % --- calculate horizontal fluxes ---
            ci = ismember(c_vic, c);
            [row, col] = find(ci);
            ci_n = c_vic(unique(row), :);
            ci_unique = unique(ci_n);
            for j=1:numel(ci_unique)
               cij_unique = ci_unique(j);
               [rj, cj] = find(ci_n == cij_unique);
               cij = ci_n(rj, :);
               vGh = 0;
               for k=1:size(cij, 1)
                  vGh = vGh + (3-2*k).*vG(cij(k, 3-cj(k))); % 3-cj(k) ensures opposite column is chosen (i.e. neighboring cell). 3-2k ensures fluxes OUT of column have consistent sign
               end
               %vG_veHorz(find(c == cij_unique)) = vGh;
               vG_veHorz(c == cij_unique) = vGh; % set flux into correct index (corresponding position for vG_veTrans)
            end
            % -----------------------------------
            
            Sn_c = value(Sn(c));
            Sn0_c = value(Sn0(c));
            sGnve0_c = value(sGnve0(c));
            
            veBottomFilled = (B == b & T ~= t) & (vG_veTrans > 0 & Sn_c < snr); % bottom of col is filled but co2 hasn't reached top yet => not in VE!
            %veTopEmpty = (B ~= b & T == t) & vG_veTrans == 0; % top of ve column empty => VE assumption not hold yet!       
            r = max(abs(vG_veHorz) ./ (abs(vG_veHorz) + vG_veTrans), 0); % max to avoid nan
            
            h = r.*(Sn0_c - sGnve0_c).*H./(1-swr);
            h_max = sGnve0_c.*H./snr + (1-r).*(Sn0_c - sGnve0_c).*H./snr + h;
            veValid = (B == b & T ~= t) & ((vG_veTrans > 0 & Sn_c >= snr & vG_veHorz == 0) | ... % fully residual filled and no horizontal fluxes from neighbors
                                            (vG_veTrans > 0 & h_max >= H*(1-1e-4) & abs(vG_veHorz) > 0 )); % height of mobile + residual part exceeds height of column (i.e. bottom residual and top mobile part "merges")
                                            
            %c_veNot = intersect(c(veBottomFilled), c(veTopEmpty)); % internal ve columns where ve assumption not satisfied
            c_veNot = c(veBottomFilled);
            c_ve = c(veValid);                                  
            
            c_veNot = c_vic(ismember(c_vic, c_veNot));           
            if ~isempty(c_veNot)
                c_VE_Not = cat(2, c_VE_Not, c_veNot);
            end
            
            c_ve = c_vic(ismember(c_vic, c_ve));
            if ~isempty(c_ve)
                c_VE = cat(2, c_VE, c_ve);                            
            end
            
            % Horizontal fluxes for NVEs
            veHorz = veBottomFilled & (abs(vG_veHorz) > 0) & ~veValid; % & Sn_c > Sn0_c;
            c_veHorz = c(veHorz);
            
            if ~isempty(c_veHorz)
                c_VE_Horz = cat(2, c_VE_Horz, c_veHorz);
                R = cat(2, R, r(veHorz));
            end
        end
    end
    
    c_VE_Not = unique(cell2mat(c_VE_Not));
    c_VE = unique(cell2mat(c_VE));
    c_VE_Horz = unique(cell2mat(c_VE_Horz));  
    R = cell2mat(R(:));
end

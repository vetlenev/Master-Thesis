function [c_VE_Not, c_VE, c_VE_Horz] = getResidualFilledCells_test(model, Sn, vG, Sn0, sGnve0)
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
    veTransition = op.connections.veToFineConn | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;
    veInternalConn = op.connections.veInternalConn;
    veHorzConn = op.connections.veTransitionHorizontalConn;
    
    c_VE_Not = {};
    c_VE = {};
    c_VE_Horz = {};
    n = op.N; 
    %n = CG.faces.neighbors;
    cn = op.N(veTransition, :);
    c_vihc = op.N(veInternalConn | veHorzConn, :);
    vG_veTrans = vG(ismember(n, cn, 'rows'));
    % NB: Get ONE horizontal flux value per cn-cell
    vG_veHorz = vG(ismember(n, c_vihc, 'rows'));
    
    for idx=1:2 % CONSIDER REMOVING LOOP: veBottomFilled and veTopEmpty will never both occur as it stand now!
        c = cn(:, idx);        
        isVE_c = ~isFine(c);
           
        if any(isVE_c) % reconstruct specifically for ve col transitioning to a lower layer
            t = CG.cells.topDepth(c);
            T = op.connections.faceTopDepth(veTransition, idx);
            b = CG.cells.bottomDepth(c);
            B = op.connections.faceBottomDepth(veTransition, idx);
            H = b - t;
            
            Sn_c = value(Sn(c));
            Sn0_c = value(Sn0(c));
            sGnve0_c = value(sGnve0(c));
            %veBottomFilled = (B == b & T ~= t) & Sn_c > 1e-4; % transition from UPPER ve column to LOWER fine cells, not from LOWER ve col to UPPER fine cells, and plume comes from lower layer.
            veBottomFilled = (B == b & T ~= t) & (vG_veTrans > 0 & Sn_c < snr); % bottom of col is filled but co2 hasn't reached top yet => not in VE!
            %veTopEmpty = (B ~= b & T == t) & vG_veTrans == 0; % top of ve column empty => VE assumption not hold yet!       
            h_max = sGnve0_c.*H./snr + (Sn0_c - sGnve0_c).*H./(1-swr);
            veValid = (B == b & T ~= t) & ((vG_veTrans > 0 & Sn_c >= snr & vG_veHorz == 0) | ... % fully residual filled and no horizontal fluxes from neighbors
                                            (vG_veTrans > 0 & h_max >= H-1e-4 & vG_veHorz > 0 )); % height of mobile + residual part exceeds height of column (i.e. bottom residual and top mobile part "merges")
                                            
            %c_veNot = intersect(c(veBottomFilled), c(veTopEmpty)); % internal ve columns where ve assumption not satisfied
            c_veNot = c(veBottomFilled);
            c_ve = c(veValid);
            
            c_vic = op.N(veInternalConn, :);            
            
            c_veNot = c_vic(ismember(c_vic, c_veNot));           
            if ~isempty(c_veNot)
                c_VE_Not = cat(2, c_VE_Not, c_veNot);
            end
            
            c_ve = c_vic(ismember(c_vic, c_ve));
            if ~isempty(c_ve)
                c_VE = cat(2, c_VE, c_ve);                            
            end
            
            % Horizontal fluxes for NVEs
            veHorz = veBottomFilled & (vG_veHorz > 0) & ~veValid;
            c_veHorz = c(veHorz);
            
            if ~isempty(c_veHorz)
                c_VE_Horz = cat(2, c_VE_Horz, c_veHorz);
            end
        end
    end
    
    c_VE_Not = unique(cell2mat(c_VE_Not));
    c_VE = unique(cell2mat(c_VE));
    c_VE_Horz = unique(cell2mat(c_VE_Horz));    
       
end

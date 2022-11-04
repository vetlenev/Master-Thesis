function [c_VE_Not, c_VE] = getResidualFilledCells(model, Sn, vG)
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
    
    c_VE_Not = {};
    c_VE = {};
    n = op.N; 
    %n = CG.faces.neighbors;
    cn = op.N(veTransition, :);
    vG_veTrans = vG(ismember(n, cn, 'rows'));
    
    for idx=1:2 % CONSIDER REMOVING LOOP: veBottomFilled and veTopEmpty will never both occur as it stand now!
        c = cn(:, idx);        
        isVE_c = ~isFine(c);
           
        if any(isVE_c) % reconstruct specifically for ve col transitioning to a lower layer
            t = CG.cells.topDepth(c);
            T = op.connections.faceTopDepth(veTransition, idx);
            b = CG.cells.bottomDepth(c);
            B = op.connections.faceBottomDepth(veTransition, idx);              
            
            Sn_c = value(Sn(c));         
            %veBottomFilled = (B == b & T ~= t) & Sn_c > 1e-4; % transition from UPPER ve column to LOWER fine cells, not from LOWER ve col to UPPER fine cells, and plume comes from lower layer.
            veBottomFilled = (B == b & T ~= t) & (vG_veTrans > 0 & Sn_c < snr); % bottom of col is filled but co2 hasn't reached top yet => not in VE!
            %veTopEmpty = (B ~= b & T == t) & vG_veTrans == 0; % top of ve column empty => VE assumption not hold yet!       
            veValid = (B == b & T ~= t) & (vG_veTrans > 0 & Sn_c >= snr); % residual plume has reached top -> from prf to frf -> VE assumption holds
                                            
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
        end
    end
    
    c_VE_Not = unique(cell2mat(c_VE_Not));
    c_VE = unique(cell2mat(c_VE));             
       
end

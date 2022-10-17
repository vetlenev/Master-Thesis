function [c_VE_Not, c_VE, snMax] = getResidualFilledCells(model, swr, snr, Sn, snMax, vG)
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
    
    tt = CG.cells.topDepth;
    %TT = CG.parent.cells.topDepth;

    bb = CG.cells.bottomDepth;
    %BB = CG.parent.cells.bottomDepth;
    H = CG.cells.height;
    
    %[a_M, a_R, sG] = getGasSatFromHeight(TT, tt, BB, bb, h, h_max, swr, snr);
    % ------------------- 
      
    % Modify for ve bottom transition
    veTransition = op.connections.veToFineConn | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;
    veInternalConn = op.connections.veInternalConn; 
     
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
            c_veNot_bool = ismember(c_vic, c_veNot);
            c_VE_Not = c_vic(c_veNot_bool);
            c_ve_bool = ismember(c_vic, c_ve);
            c_VE = c_vic(c_ve_bool);
            
%             for i=1:2
%                 c_vic_i = c_vic(:,i);
%                 c_veNot_i = c_veNot_bool(:,i);
%                 c_ve_i = c_ve_bool(:,i);
%                 c_VE_Not{i} = c_vic_i(c_veNot_i); % get ve cells whose bottom virtual cell satisfies sG > 1e-4
%                 c_VE{i} = c_vic_i(c_ve_i);
%             end                   
        end
    end
    
    c_VE_Not = unique(c_VE_Not);
    c_VE = unique(c_VE);
    c_vic = unique(c_vic); % all unique internal ve cells
    Sn_vic = Sn(c_vic);% only update internal ve columns
    snMax_vic = snMax(c_vic);
    
    if ~isempty(c_VE) % once residual plume reaches top, set sgMax to maximum and use standard formulas       
        snMax(c_VE) = 1-swr; % would otherwise be discontinuous transition from H to snr*H
    end  
       
    if ~isempty(c_VE_Not) % only update reconstructed saturation if co2 originates from bottom of any ve column
        % --- Partly residual filled (frf) ---        
        mask_VE_Not = ismember(c_vic, c_VE_Not);       
        % choose max of current and previous max
        Sn_vic(mask_VE_Not) = max(Sn_vic(mask_VE_Not), snMax_vic(mask_VE_Not));
        % update global coarse saturation
        %sG(c_VE_Not) = Sn_vic(mask_VE_Not);       
    end                  
       
end

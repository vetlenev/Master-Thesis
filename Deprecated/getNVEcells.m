function [c_VE_Not, c_VE] = getNVEcells(model, sG, vG)
    % Get NVE cells and newly transformed NVE->VE cells,
    % based on coarse saturation and bottom flux.
    CG = model.G;
    op = model.operators;
    % Modify for ve bottom transition
    veTransition = op.connections.veToFineConn | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;
    veInternalConn = op.connections.veInternalConn; 
      
    %n = CG.faces.neighbors;
    n = op.N;
    cn = op.N(veTransition, :);
    vG_veTrans = vG(ismember(n, cn, 'rows'));
    
    %vG_any = max(vG_veTrans > 0, vG_any); % if any vG goes from zero to non-zero, vG_any updated to 1
    
    for idx=1:2 % loop over neighbor set
        c = cn(:, idx);         
        isVE_c = ~isFine(c);
        
        if any(isVE_c) % reconstruct specifically for ve col transitioning to a lower layer
            t = CG.cells.topDepth(c);
            T = op.connections.faceTopDepth(veTransition, idx);
            b = CG.cells.bottomDepth(c);
            B = op.connections.faceBottomDepth(veTransition, idx);                                                             
                                    
            Sn_c = value(sG(c));
            veBottomConnIdx = (B == b & T ~= t) & (vG_veTrans > 0 & Sn_c < snr); % transition from UPPER ve column to LOWER fine cells, not from LOWER ve col to UPPER fine cells, and plume comes from lower layer.
            %veBottomConnIdx = (B == b & T ~= t) & (abs(vG_veTrans) > 0 | vG_any); % nonzero flux of co2 from bottom
            veValid = (B == b & T ~= t) & (vG_veTrans > 0 & Sn_c >= snr);
            
            c_veNot = c(veBottomConnIdx);
            c_ve = c(veValid);
   
            c_vic = op.N(veInternalConn, :);
            c_veNot_bool = ismember(c_vic, c_veNot);
            c_ve_bool = ismember(c_vic, c_ve);
            c_VE_Not = c_vic(c_veNot_bool);
            c_VE = c_vic(c_ve_bool);
                 
        end
    end
    
    c_VE_Not = unique(c_VE_Not);
    c_VE = unique(c_VE);   
end
function [vW, vG, mobW, mobG, upcw, upcg, h] = computeHybridFluxesVEpc(model, pW, sG, muW, muG, rhoW, rhoG, trans, varargin)
% Internal function - computes interior fluxes for the hybrid VE models

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    op = model.operators;
    G = model.G;
    g = norm(model.gravity);
    
    H = G.cells.height;    
    h = H.*sG;
    isFine = G.cells.discretization  == 1;    
    
    [pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, sG, h, H, rhoW, rhoG, muW, muG, isFine, isFine);   
    
    if isa(model, 'ThreePhaseCompositionalModel')
        sW = 1-sG;
        rhoWf = op.faceAvg(rhoW.*sW)./max(op.faceAvg(sW), 1e-8);
        rhoGf = op.faceAvg(rhoG.*sG)./max(op.faceAvg(sG), 1e-8);
    else
        rhoWf = op.faceAvg(rhoW);
        rhoGf = op.faceAvg(rhoG);
    end

    z = (G.cells.bottomDepth + G.cells.topDepth)/2;
    gdz =  g.*op.Grad(z);
    
    vIc = op.connections.veInternalConn;
    n1 = op.N(vIc, 1);
    n2 = op.N(vIc, 2);
    
    
    [gdz_g, gdz_w] = deal(gdz);
    gdz_g(vIc) = g*(G.cells.topDepth(n2) - G.cells.topDepth(n1));
    gdz_w(vIc) = g*(G.cells.bottomDepth(n2) - G.cells.bottomDepth(n1));        
    
    dpG   = op.Grad(pG) - rhoGf .* gdz_g;
    upcg  = (value(dpG)<=0);
    vG = -op.faceUpstr(upcg, mobG) .* trans .* dpG;

    dpW   = op.Grad(pW) - rhoWf .* gdz_w;
    upcw  = (value(dpW)<=0);
    vW = -op.faceUpstr(upcw, mobW) .* trans .* dpW;

    % Treat transition between different discretizations
    transition_ve = model.operators.connections.veTransitionHorizontalConn & ~model.operators.connections.veToFineConn;

    if any(transition_ve)
        [vW(transition_ve), vG(transition_ve),...
         upcw(transition_ve), upcg(transition_ve)] = computeTransitionFluxVE(model, pW, pG, h, rhoW, rhoG, muW, muG, transition_ve, true);
    end
    % Treat transition between ve zones in vertical direction (for diffuse
    % leakage)
    transition_vertical = model.operators.connections.veToFineConn | ...
                          model.operators.connections.veTransitionVerticalConn & op.T > 0;
  
    if any(transition_vertical)        
        [vW(transition_vertical), vG(transition_vertical),...
         upcw(transition_vertical), upcg(transition_vertical)] = computeTransitionFluxVE(model, pW, pG, h, rhoW, rhoG, muW, muG, transition_vertical, false);
    end   
    
end

function [pW, sG, h, H, rhow, rhog, muW, muG, isFine] = getTransitionValuesVE_coarse(model, pW, h_global, index, subs, rhoW, rhoG, muW, muG)    
    c = model.operators.N(subs, index);
    
    % Initialize variables for full domain
%     t = UtilFunctions.initNanADI(h_global); T = UtilFunctions.initNanADI(h_global);
%     b = UtilFunctions.initNanADI(h_global); B = UtilFunctions.initNanADI(h_global);
%     H = UtilFunctions.initNanADI(h_global);
    rhow = UtilFunctions.initNanADI(h_global); rhog = UtilFunctions.initNanADI(h_global);
    %muw = UtilFunctions.initNanADI(h_global); mug = UtilFunctions.initNanADI(h_global);          
    pW_full = UtilFunctions.initNanADI(h_global);
    
    isFine = model.G.cells.discretization == 1;
    
    % Fill values only for transition cells
    t = model.G.cells.topDepth(c);
    %T(c) = model.operators.connections.faceTopDepth(subs, index);  
    T = model.operators.connections.faceTopDepth(subs, index);
    b = model.G.cells.bottomDepth(c);
    B = model.operators.connections.faceBottomDepth(subs, index);
    H = model.operators.connections.faceHeight;
    %H = model.operators.connections.faceHeight(subs, index);

    sG = getGasSaturationFromHeight(T, t, B, b, h_global, model, c);
    %sG = getGasSaturationFromHeight(T(c), t(c), B(c), b(c), h_global(c), model, c);
    h = sG(c).*H; 
    %h = sG.*H;
    
    g = norm(model.gravity);
    rhow(c) = rhoW(c);
    rhog(c) = rhoG(c);
    
    %muw(c) = muW(c);
    %mug(c) = muG(c);
    
    %pW_full(c) = pW(c);
%     pW = getWaterPressureFromHeight(B(c), t(c), B(c), b(c), h_global(c), ...
%                                     pW_full(c), g, rhow(c), rhog(c), model, c);
    pW = getWaterPressureFromHeight(B, t, B, b, h_global, ...
                                pW(c), g, rhow(c), rhog(c), model, c);
end

function [pW, sG, h, H, rhow, rhog, muW, muG, isFine] = getTransitionValuesVE_fine(model, pW, h_global, index, subs, rhoW, rhoG, muW, muG)      
    c = model.operators.N(subs, index);
    % Initialize variables for full domain
%     t = UtilFunctions.initNanADI(h_global); T = UtilFunctions.initNanADI(h_global);
%     b = UtilFunctions.initNanADI(h_global); B = UtilFunctions.initNanADI(h_global);
%     H = UtilFunctions.initNanADI(h_global); C = UtilFunctions.initNanADI(h_global);
%     rhow = UtilFunctions.initNanADI(h_global); rhog = UtilFunctions.initNanADI(h_global);
%     muw = UtilFunctions.initNanADI(h_global); mug = UtilFunctions.initNanADI(h_global);          
    pW_f = UtilFunctions.initNanADI(h_global);
    h = UtilFunctions.initNanADI(h_global);
    
    %isFine = model.G.cells.discretization(c) == 1;
    isFine = model.G.cells.discretization == 1;
       
    % Fill values only for transition cells    
    t = model.G.cells.topDepth(c);
    T = model.operators.connections.faceTopDepth(subs, index);  
    b = model.G.cells.bottomDepth(c);
    B = model.operators.connections.faceBottomDepth(subs, index);
    H = model.operators.connections.faceHeight(subs, height);
    
    %sG = getGasSaturationFromHeight(T, t, B, b, h_global(c));
    %sG = getGasSaturationFromHeight(T(c), t(c), B(c), b(c), h_global, model, c);
    sG = getGasSaturationFromHeight(T, t, B, b, h_global, model, c);
    h(c) = sG(c).*H;
    
    g = norm(model.gravity);
    %rhow(c) = rhoW(c);
    %rhog(c) = rhoG(c);
    C = (T + B)/2;
    pW_f(c) = pW(c);
    %muw(c) = muW(c);
    %mug(c) = muG(c);
   
    %pW_c = getWaterPressureFromHeight(C, t, B, b, h_global(c), pW_f, g, rhow, rhog);
%     pW_c = getWaterPressureFromHeight(C(c), t(c), B(c), b(c), h_global, ...
%                                         pW_f(c), g, rhoW(c), rhoG(c), model, c);
    pW_c = getWaterPressureFromHeight(C, t, B, b, h_global, ...
                                        pW_f(c), g, rhoW(c), rhoG(c), model, c);

    pW = isFine.*pW_f + ~isFine.*pW_c;
end


function [vW, vG, upcw, upcg] = computeTransitionFluxVE(model, pW, pG, h, rhoW, rhoG, muW, muG, vtc, treatAsCoarse, varargin)
    op = model.operators;
    G = model.G;
    trans = model.operators.T(vtc);
    
    % Capillary pressure for all cells
    pcWG = pG - pW;
    
    % Consistent with rest of MRST...
    grad = @(l, r) r - l;
    upstr = @(flag, l, r) flag.*l + ~flag.*r;
    rhoWf = op.faceAvg(rhoW);
    rhoGf = op.faceAvg(rhoG);
    
    rhoWf = rhoWf(vtc);
    rhoGf = rhoGf(vtc);
    g = norm(model.gravity);
    t1 = op.connections.faceTopDepth(vtc, 1);
    t2 = op.connections.faceTopDepth(vtc, 2);

    b1 = op.connections.faceBottomDepth(vtc, 1);
    b2 = op.connections.faceBottomDepth(vtc, 2);
    isFine = G.cells.discretization == 1;
    if treatAsCoarse
        % Treat fine cells as VE cells
        if any(isFine)
            % Fine cells: Calculate the water pressure at the bottom of the
            % cell by assuming hydrostatic equilibrium within each fine cell.
            % This is then backwards corrected in the "VE" way by further calls
            % to the compute routines.

            T = G.cells.topDepth(isFine);
            B = G.cells.bottomDepth(isFine);
            H = G.cells.height(isFine);
            
            cz = (T + B)/2;
            dz = B - cz;
            dpFine = g.*rhoW(isFine).*dz;
            pW(isFine) = pW(isFine) - dpFine;
        end    

        [pW_l, sG_l, h_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesVE_coarse(model, pW, h, 1, vtc, rhoW, rhoG, muW, muG);
        [pW_r, sG_r, h_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesVE_coarse(model, pW, h, 2, vtc, rhoW, rhoG, muW, muG);
       
        isFine_l = false;
        isFine_r = false;
        gdz_g = g*grad(t1, t2);
        gdz_w = g*grad(b1, b2);
    else
        % Treat coarse cells as fine cells              
        [pW_l, sG_l, h_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesVE_fine(model, pW, h, 1, vtc, rhoW, rhoG, muW, muG);       
        [pW_r, sG_r, h_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesVE_fine(model, pW, h, 2, vtc, rhoW, rhoG, muW, muG);

        isFine_l = true;
        isFine_r = true;
        c1 = (t1 + b1)/2;
        c2 = (t2 + b2)/2;
        gdz_g = g*grad(c1, c2);
        gdz_w = gdz_g;
    end
    
    nn = model.operators.N(vtc, :);
    n_l = nn(:,1); n_r = nn(:,2);
    %nn_fine = model.G.cells.discretization(nn) == 1;              

    [pW_l, pG_l, mobW_l, mobG_l] = evaluatePropertiesVE(model, pW_l, sG_l, h_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l, isFine_l, n_l, pcWG);
    [pW_r, pG_r, mobW_r, mobG_r] = evaluatePropertiesVE(model, pW_r, sG_r, h_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r, isFine_r, n_r, pcWG);            
    
    pW_l = pW_l(n_l); pG_l = pG_l(n_l); 
    mobW_l = mobW_l(n_l); mobG_l = mobG_l(n_l);
    pW_r = pW_r(n_r); pG_r = pG_r(n_r); 
    mobW_r = mobW_r(n_r); mobG_r = mobG_r(n_r);        
    
    if ~isempty(varargin)
        transFunc = varargin{1};
        pE = varargin{2};
        pcWG_l = pG_l - pW_l;
        pcWG_r = pG_r - pW_r;
        % -------- Or choose to be cell with max sG
        %pcWG = max(pcWG_l, pcWG_r); % max pc is the cell below the interface ???
        if mean(sG_l(n_l)) > mean(sG_r(n_r))
           pcWG = pcWG_l;
        else
           pcWG = pcWG_r;
        end
        % --------
        trans = transFunc(trans, pcWG, pE); % transmissibility as function of capillary pressure
    end
    
    dpG   = grad(pG_l, pG_r) - rhoGf .* gdz_g;
    upcg  = (value(dpG)<=0);   
    vG = -upstr(upcg, mobG_l, mobG_r) .* trans .* dpG;

    dpW   = grad(pW_l, pW_r) - rhoWf .* gdz_w;
    upcw  = (value(dpW)<=0);
    vW = -upstr(upcw, mobW_l, mobW_r) .* trans .* dpW;
end


function [pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, sG, h, H, rhoW, rhoG, muW, muG, isFine, nn_fine, varargin)
    g = norm(model.gravity);
    f = model.fluid;
    sW = 1 - sG;
    isVE = ~isFine;
    
    if ~isempty(varargin)
       pcWG_all = varargin{1}; % store computed values for entire domain
    end
    
    if isfield(f, 'pcWG')        
        pcWG = f.pcWG(sG, nn_fine);
    else
        pcWG = 0;
    end
    
    pcWG = pcWG + isVE.*(h.*(rhoW - rhoG) - H.*rhoW).*g;
    pG = pW + pcWG;
    
    if ~isempty(varargin) % retrieve non-transition values
       indexMap = (1:model.G.cells.num)';
       non_transition = setdiff(indexMap, nn_fine);
       transition = setdiff(indexMap, non_transition);
       pcWG(non_transition) = pcWG_all(non_transition); % extract only non-transition values
       % Gas pressure
       pG(non_transition) = pW(non_transition) + pcWG(non_transition);   
    end
    
    % Mobility
    if isprop(model, 'EOSModel')
        krw = f.krO(sW);
    else
        krw = f.krW(sW);
    end
    mobW = (isVE.*sW + isFine.*krw)./muW;
    mobG = (isVE.*sG + isFine.*f.krG(sG))./muG;
end

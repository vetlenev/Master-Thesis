function [masses, masses_0] = massTrappingFine_Other(Gh, cFine, p, sG, sW, rock, fluidADI, varargin)
% Compute the trapping status distribution of CO2 in FINE cells not part of
% a subgrid for a sealing layer.
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gh         - hybrid grid
%   cFine      - remaining fine cells
%   p          - pressure, one value per cell of grid
%   sW         - water saturation, one value per cell of HYBRID grid
%   sG         - gas saturation, one value per cell of grid
%   rock       - rock parameters corresponding to 'Gt'
%   fluidADI   - ADI fluid object (used to get densities and compressibilities)
%   varargin   - optional parameters/value pairs.  This currently only
%                includes the option 'rs', which specifies the amount of
%                dissolved CO2 (in its absence, dissolution is ignored).
%
% RETURNS:
%   masses - vector with 7 components, representing:
%            masses[1] : mass of dissolved gas, per cell
%            masses[2] : mass of gas that is both structurally and residually trapped
%            masses[3] : mass of gas that is residually (but not structurally) trapped
%            masses[4] : mass of non-trapped gas that will be residually trapped
%            masses[5] : mass of structurally trapped gas, not counting the gas that 
%                        will eventually be residually trapped
%            masses[6] : mass of subscale trapped gas (if 'dh' is nonempty)
%            masses[7] : mass of 'free' gas (i.e. not trapped in any way)
%   masses_0 (optional) - masses in terms of one value per grid cell of Gt

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

    opt.rs = 0;
    opt    = merge_options(opt, varargin{:});       
    
    % Extracting relevant information from 'sol'
    sw=fluidADI.krPts.w(1);%liquid residual saturation (scalar)
    sr=fluidADI.krPts.g(1);%gas residual saturation (scalar)
    SW = sW(cFine);
    SG = sG(cFine);
    rs = opt.rs;    
    p = p(cFine);
    
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(p);
    end
    pv       = rock.poro(cFine) .* Gh.cells.volumes(cFine) .* pvMult; % think Gt.cells.volumes is fine -> it's just the surface area of each top surface cell
    if isfield(rock,'ntg')
       pv = rock.poro(cFine) .* Gh.cells.volumes(cFine) .* rock.ntg(cFine) .* pvMult;
       %effective area accounting for possible net-to-gross data
    end
    rhoCO2   = fluidADI.rhoGS .* fluidADI.bG(p);
    gasPhase = sum(pv .* (rhoCO2 .* SG));                   
                      
    % No structural trapping for remaining cells!
    resStruc  = 0;      % trapped, res
    freeStruc = 0;      % trapped, non-res 
        
    plumeSG = SG(SG > sr); % part of mobile plume
    freeRes   = sum(pv(SG > sr).*rhoCO2(SG > sr).*sr);                                     % non-trapped, flowing, rrrrres
    freeMov   = sum(pv(SG > sr).*rhoCO2(SG > sr).*max(plumeSG - sr, 0));                          % non-trapped, flowing, non-res 
    resSG = SG(SG <= sr); % part of immobile plume       
    resTrap = sum(pv(SG <= sr) .* rhoCO2(SG <= sr) .* resSG);
    
    resDis    = fluidADI.rhoGS .* sum(pv.* (rs .* fluidADI.bW(p) .* SW)); % dissolved
    subtrap   = 0;
          
    masses    = max([value(resDis), value(resStruc), value(resTrap), ...
                     value(freeRes), value(freeStruc), value(subtrap), ...
                     value(freeMov)], 0); % may be ADI variables
         
    if(abs(sum(masses(2:end))-gasPhase) > 1e-3 * gasPhase)
        %disp('There is a mismatch between mass calculations');
        %fprintf('Mass: %d. Gas: %d\n', abs(sum(masses(2:end))-gasPhase), 1e-3 * gasPhase);       
    end
    
    if nargout > 1
        % values one per cell of grid Gt
        plumeSG_0 = zeros(numel(SG), 1);
        resSG_0 = zeros(numel(SG), 1);
        plumeSG_0(SG > sr) = SG(SG > sr);
        resSG_0(SG <= sr) = SG(SG <= sr);
             
        resStruc_0  = 0;
        freeStruc_0 = 0;
        freeRes_0   = pv(SG > sr).*rhoCO2.*sr;
        freeMov_0   = pv(SG > sr).*rhoCO2.*max(plumeSG_0 - sr, 0);
        resTrap_0   = pv(SG <= sr) .* rhoCO2 .* resSG_0;
        
        resDis_0    = fluidADI.rhoGS .* pv .* (rs .* fluidADI.bW(p) .* SW);
        subtrap_0   = 0;

        masses_0 = [{resDis_0}, {resStruc_0}, {resTrap_0}, {freeRes_0}, ...
                    {freeStruc_0}, {subtrap_0}, {freeMov_0}]; % may be ADI vars
    end
end

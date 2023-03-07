function f = assignSGFN(f, sgfn, reg)
    %
    % LS: Modified entirely based on assignSGOF. This avoids the usage of
    % interpReg and getRegMap (deprecated) and is compatible with
    % assignEHYSTR.
    %
%     if reg.sat == 1                     % single region (drainage)
%         regs = 1;                             
%         imb   = false;
%     elseif reg.sat > numel(reg.SATINX)  % both dr and imb regions
%         regs = reg.sat/2;  
%         regi = reg.sat/2;
%         imb   = true;
%     else                                % multiple dr regions
%         regs = reg.sat;  
%         imb   = false;
%     end
    
    if reg.sat > numel(reg.SATINX)
        imb   = true;
    else
        imb = false;
    end
    regs = reg.sat;
    [f.krG, f.pcOG] = deal(cell(1, regs));
    for i = 1:regs
        SGFN      = extendTab(sgfn{i});
        SG        = SGFN(:, 1);
        f.krG{i}  = @(sg) interpTable(SG, SGFN(:, 2), sg);
        f.pcOG{i} = @(sg) interpTable(SG, SGFN(:, 3), sg);
    end
    
    if imb
        nregi = regs/2;
        for i = 1:nregi
            SGFN        = extendTab(sgfn{i+nregi});
            SG          = SGFN(:, 1);
            spos        =  SG > 0;
            kr0         = SGFN(:,2) == 0;
            f.sgtmax(i) = SGFN(all([spos kr0], 2), 1)';
        end
        assert(numel(f.sgtmax) == nregi)
    end
end

% function f = assignSGFN(f, sgfn, reg)
%    cfun = @(f) cellfun(f, sgfn, 'UniformOutput', false);
% 
%    % Compute tables (static data)
%    Tkrg  = extendTab( cfun(@(x) x(:, [1, 2])) );
%    Tpcog = extendTab( cfun(@(x) x(:, [1, 3])) );
% 
%    % Region mapping
%    regmap = @(sg, varargin) ...
%       getRegMap(sg, reg.SATNUM, reg.SATINX, varargin{:});
% 
%    % Region interpolator
%    ireg = @(T, sg, varargin) interpReg(T, sg, regmap(sg, varargin{:}));
% 
%    f.krG  = @(sg, varargin) ireg(Tkrg , sg, varargin{:});
%    f.pcOG = @(sg, varargin) ireg(Tpcog, sg, varargin{:});
%    
%    % LS: Added to assign bounding imbibition curve to f if there is one.
%    if reg.sat > numel(reg.SATINX)
%        if reg.sat == 2
%            reg.IMBNUM = 2;
%        end       
%        regmapi = @(sg, varargin) ...
%            getRegMap(sg, reg.IMBNUM, reg.IMBINX, varargin{:});
%        
%        iregi = @(T, sg, varargin) interpReg(T, sg, regmapi(sg, ...
%                                                            varargin{:}));
%        
%        % Bounding imbibition and Pc curves                                                
%        f.krGib  = @(sg, varargin) iregi(Tkrg , sg, varargin{:});
%        f.pcOGib = @(sg, varargin) iregi(Tpcog, sg, varargin{:});
%        
%        % Maximum trapped gas saturation
%        Tkrgi = cell2mat(Tkrg(unique(reg.IMBNUM)));  
%        spos  =  Tkrgi(:,1) > 0;     
%        kr0   = Tkrgi(:,2) == 0;
%        f.sgtmax = Tkrgi(all([spos kr0], 2), 1)';
%        assert(numel(f.sgtmax) == reg.sat/2)
%    end % LS <
%    
% end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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

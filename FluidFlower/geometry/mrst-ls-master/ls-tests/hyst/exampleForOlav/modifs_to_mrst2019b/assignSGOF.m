function f = assignSGOF(f, sgof, reg)
    % Modified by LS to account for nonwetting phase relperm hysteresis.
    [f.krG, f.krOG, pcOG, f.krPts.g, f.krPts.og, hasPC] = getFunctions(f, sgof, reg);
    if hasPC
        f.pcOG = pcOG;
    end
%     if imb
%         f.krGib    = krGib;
%         f.krOGib   = krOGib;
%         f.krPts.gi = ptsi;
%         if hasPCi
%             f.pcOGib = pcOGib;
%         end
%     end
end

function [krG, krOG, pcOG, pts, pts_o, hasPC] = getFunctions(f, SGOF, reg)

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
    regs = reg.sat;
    [krG, krOG, pcOG] = deal(cell(1, regs));
    [pts, pts_o] = deal(zeros(regs, 4));
    hasPC = false;
    for i = 1:regs
        if isfield(f, 'krPts')
            swcon = f.krPts.w(i, 1);
        else
            warning('No relperm points found in assignment of SGOF.');
            swcon = 0;
        end
        [pts(i, :), pts_o(i, :)] = getPoints(SGOF{i}, swcon);
        sgof = extendTab(SGOF{i});
        SG = sgof(:, 1);
        pc = sgof(:, 4);
        hasPC = hasPC || any(pc ~= 0);
        krG{i} = @(sg) interpTable(SG, sgof(:, 2), sg);
        krOG{i} = @(so) interpTable(SG, sgof(:, 3), 1-so-swcon);
        pcOG{i} = @(sg) interpTable(SG, pc, sg);
    end
    
%     hasPCi                  = false;
%     [krGib, krOGib, pcOGib] = deal(cell(1, regi));
%     ptsi                    = deal(zeros(regs, 4));
%     if imb
%         for i = 1:regi
%             ptsi(i, :)  = getPoints(SGOF{i+regs}, swcon);
%             sgofi       = extendTab(SGOF{i+regs});
%             SG          = sgofi(:, 1);
%             pci         = sgofi(:, 4);
%             hasPCi      = hasPCi || any(pci ~= 0);
%             krGib{i}  = @(sg) interpTable(SG, sgofi(:, 2), sg);
%             krOGib{i} = @(sg) interpTable(SG, sgofi(:, 3), sg);
%             pcOGib{i} = @(sg) interpTable(SG, pci, sg);
%         end
%     end
end

function [pts, pts_o] = getPoints(sgof, swcon)
    % Connate gas saturation
    pts = zeros(1, 4);
    pts(1) = sgof(1, 1);
    % Last mobile gas saturation
    ii = find(sgof(:,2)==0, 1, 'last');
    pts(2) = sgof(ii,1);
    % Last point
    pts(3) = sgof(end,1);
    % Maximum relperm
    pts(4) = sgof(end,2);
    
    % Get OG-scaling
    pts_o = zeros(1, 4);
    pts_o(3) = 1;
    ii = find(sgof(:,3) == 0, 1, 'first');
    pts_o(2) = 1 - sgof(ii,1) - swcon;
    % Maximum oil relperm
    pts_o(4) = sgof(1,3);
end

function v = selectSubset(v, varargin)
if ~isempty(varargin)
    v = v(varargin{2});
end
end
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

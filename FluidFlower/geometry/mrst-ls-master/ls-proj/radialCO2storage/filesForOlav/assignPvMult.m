function f = assignPvMult(f, cP, idreg)
% This function takes an existing fluid object and assigns a pore volume
% multiplier to it taking into account the inital pressure (and not the
% reference pressure indicated for the cr in the .DATA input file)
%
% INPUTS:
% f:     fluid object
% cP:    pore (rock) compressibility in 1/Pa
% idreg: vector array of dimensions nx1 where n is the total number of grid
%        cells. Each entry is a single double value indicating which region
%        that cell pertains to.
%
% OUTPUTS:
% f:    fluid object with f.pvMultR assigned
%
%

% check if pvMultR is already an existing field and remove it
if isfield(f, 'pvMultR')
    f = rmfield(f, 'pvMultR');
    warning('f.pvMultR was already existent. Now overwritten.')
end

nreg = max(idreg);
if nreg == 1
    f.pvMultR = @(p, p0) pvMult(p, cP, p0);
else
    f.pvMultR = cell(1, nreg);
    for i = 1:nreg
        f.pvMultR{i} = @(p, p0) pvMult(p, cP(i), p0);
    end
end

% function v = pvMult(p, cR, pRef)  % Linearized (simplified model)
%     v = 1 + cR.*(p-pRef);
% end

function v = pvMult(p, cP, p0)      % actual integration of  (1/phi)*(dphi/dp)
    v = 1.*exp(cP.*(p-p0));
end


end

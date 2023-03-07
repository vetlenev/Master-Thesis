function [fdip, fT, phi_c, phi_s, SSFc, len, permS, permC, nSim] = smearParamVals(changeClay, fdip, fT, phi_c, phi_s, SSFc, len, permS, permC)
%
% repeat variations in parameters to consider all possible combinations
% within a single loop.
%
if changeClay == 0
    % permS and permC must be passed as (nSim(permS or permC), 2)
    nVals = [size(permS, 1), numel(permC, 1), numel(len), numel(SSFc), numel(phi_c), numel(phi_s), numel(fT), numel(fdip)];
    len = reshape(len, numel(len), 1);
    SSFc = reshape(SSFc, numel(SSFc), 1);
    phi_c = reshape(phi_c, numel(phi_c), 1);
    phi_s = reshape(phi_s, numel(phi_s), 1);
    fT = reshape(fT, numel(fT), 1);
    fdip = reshape(fdip, numel(fdip), 1);
elseif changeClay == 1 
    % All clay inputs (perm, phic, SSFc, maxLen) variable for each clay
    % layer.
    % Must input all of these clay inputs for each clay layer.
    % Must be passed with size (nClays, 2) or (nClays x nSim(each), 2) for permS 
    % and permC, and (1, nClays) or (nSim(each), nClays) for all others.
    nVals = [size(permS, 1)/size(SSFc, 2), size(permC, 1)/size(SSFc, 2), ...
             size(len, 1), size(SSFc, 1), size(phi_c, 1), size(phi_s, 1), ...
             size(fT, 1), size(fdip, 1)]; 
         if any(nVals < 1)
             nVals(nVals < 1) = 1;
         end
end

nSim  = prod(nVals); 
permS = repmat(permS, nSim/nVals(1), 1);
permC = repmat(repmat(permC, nVals(1), 1), nSim/prod(nVals(1:2)), 1);
len   = repmat(repmat(len, prod(nVals(1:2)), 1), nSim/prod(nVals(1:3)), 1);
SSFc  = repmat(repmat(SSFc, prod(nVals(1:3)), 1), nSim/prod(nVals(1:4)), 1);
phi_c = repmat(repmat(phi_c, prod(nVals(1:4)), 1), nSim/prod(nVals(1:5)), 1);
phi_s = repmat(repmat(phi_s, prod(nVals(1:5)), 1), nSim/prod(nVals(1:6)), 1);
fT    = repmat(repmat(fT, prod(nVals(1:6)), 1), nSim/prod(nVals(1:7)), 1);
fdip  = repmat(fdip, prod(nVals(1:7)), 1);

end
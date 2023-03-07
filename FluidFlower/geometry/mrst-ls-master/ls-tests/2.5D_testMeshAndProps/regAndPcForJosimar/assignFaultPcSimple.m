function [f] = assignFaultPcSimple(f, fault, refPc)
%
% See Zulqarnain et al., IJGGC (2018; equation 3)
%
% NOTE that during simulateScheduleAd, within (My)BlackOilPc, sG for all
% fault cells (used here) is passed in ascending order of global fault cell
% id. Hence, fperm and fporo must also be passed in ascending order of
% global fault cell id, otherwise the evaluation would be wrongly done.
%

fperm = fault.perm(:, end);                                                % maximum perm will usually be in the along (updip) direction
fporo = fault.poro;
f.pcOG{fault.satRegNum} = @(sG)   pcFault(sG, fperm, fporo, refPc);

function v = pcFault(sG, perm, poro, refPc)
    v = refPc.val(sG).*sqrt((refPc.perm*poro)./(refPc.poro*perm));
end

end
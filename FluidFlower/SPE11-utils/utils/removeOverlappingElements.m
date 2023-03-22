function array_new = removeOverlappingElements(array)
%REMOVEOVERLAPPINGELEMENTS Summary of this function goes here
%   Detailed explanation goes here
    array_new = sort(array);
    if true
        idx_nonunique = diff(array_new) == 0;
        array_remove = array_new(idx_nonunique);
        array_new(ismember(array_new, array_remove)) = [];        
    else        
        idx = [find(array_new(1:end-1) ~= array_new(2:end)); length(array_new)];
        occurences = diff([0; idx]);
        array_rem = array(idx); % remove duplicates
        new_array = array_rem(occurences == 1);
    end
end


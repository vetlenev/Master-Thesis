function [c_ix, ixs] = geo_curve_loop_read(str)

    % get curve loop index
    c_ix = regexp(str, 'Curve Loop\(([0-9]+)\)', 'tokens');
    c_ix = str2num(c_ix{1}{1});

    % get indices
    tmp = regexp(str, '{(.*)}', 'tokens');
    ixs = cellfun(@str2num, split(tmp{1}, ','));
    
end

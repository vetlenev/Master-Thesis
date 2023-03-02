function [l_ix, ixs] = geo_line_read(str)

    % get line index
    l_ix = regexp(str, 'Line\(([0-9]+)\)', 'tokens');
    l_ix = str2num(l_ix{1}{1});

    % get line endpoint indices
    tmp = regexp(str, '{(.*)}', 'tokens');
    ixs = cellfun(@str2num, split(tmp{1}, ','));
    
end

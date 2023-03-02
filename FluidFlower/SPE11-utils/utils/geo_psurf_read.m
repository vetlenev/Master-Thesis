function [name, surf_ix, ixs] = geo_psurf_read(str)

    % get name and surface index
    tmp = regexp(str, 'Physical Surface\("(.*)", ([0-9]+)\)', 'tokens');
    name = tmp{1}{1};
    surf_ix = str2num(tmp{1}{2});
    
    tmp = regexp(str, '{(.*)}', 'tokens');
    ixs = cellfun(@str2num, split(tmp{1}, ','));
    
    
end
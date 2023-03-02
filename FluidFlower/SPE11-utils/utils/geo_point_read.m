function [p_ix, coords] = geo_point_read(str)

% get point index
p_ix = regexp(str, 'Point\(([0-9]+)\)', 'tokens');
p_ix = str2num(p_ix{1}{1});

% get point coords
tmp = regexp(str, '{(.*)}', 'tokens');
coords = split(tmp{1}, ',');
coords = cellfun(@str2num, coords, 'uniformoutput', false);
coords = coords(cellfun(@isscalar, coords)); % keep only those that are numbers
coords = horzcat(coords{:});

end

function [pts, loops, facies] = parse_spe11_geo(filename)

    % open file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Unable to open file')
    end
    
    pts = [];
    lines = [];
    ltags = []; loop_entries = {};
    
    fnames = {}; ftags = []; facies = {};
    
    tline = fgetl(fid);
    while ischar(tline)
        
        keyword = regexp(tline, '^[a-zA-Z ]*', 'match');
        if ~isempty(keyword)
            keyword = keyword{1};
            
            switch keyword
              case 'Point'
                [p_ix, coords] = geo_point_read(tline);
                pts = [pts; [p_ix, coords]];
              case 'Line'
                [l_ix, ixs] = geo_line_read(tline);
                lines = [lines; [l_ix, ixs(:)']];
              case 'Curve Loop'
                [c_ix, ixs] = geo_curve_loop_read(tline);
                
                ltags = [ltags; c_ix];
                loop_entries = [loop_entries, {ixs}];
              case 'Physical Surface'
                [name, surf_ix, ixs] = geo_psurf_read(tline);
                fnames = [fnames, {name}];
                ftags = [ftags; surf_ix];
                facies = [facies; {ixs}];
              otherwise
                % do not process this line
                fprintf('Skipped reading line with keyword: "%s".\n', keyword);
            end
            
            %disp(keyword);
        end
        
        tline=fgetl(fid);
    end
    fclose(fid);    
    
    % sort and check all read tables
    pts = process_and_check(pts);
    lines = process_and_check(lines);
    [ltags, order] = process_and_check(ltags); loop_entries = loop_entries(order);
    [ftags, order] = process_and_check(ftags); fnames = fnames(order); facies = facies(order);
    
    % remove index column for points and lines
    pts = pts(:, 2:end);
    lines = lines(:, 2:end);
    
    % construct  loops
    loops = construct_loops(lines, loop_entries);
    
end


function [m, order] = process_and_check(m)
    [m, order] = sortrows(m, 1);
    assert(all(m(:,1) == [1:size(m, 1)]')); % check that all entries are ordered
                                            % consecutively from 1 upwards
    
end

function loops = construct_loops(segments, loop_entries)

    loops = {};
    for l = loop_entries(:)'
        lix = [];
        l = l{:};
        for i = l(:)'
            seg = segments(abs(i),:);
            if i > 0
                lix = [lix, seg(2)];
            else
                lix = [lix, seg(1)];
            end
        end
        loops = [loops, {lix}];
    end
end
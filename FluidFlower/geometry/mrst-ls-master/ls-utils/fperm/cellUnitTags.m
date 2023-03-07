function [G] = cellUnitTags(G, ucids, faultDiv)
%
%
%

G.cells.unit = nan(G.cells.num, 1);
G.cells.Vcl  = nan(G.cells.num, 1);
if faultDiv == 1                                                            % separate ids for fault segments juxtaposing different units
    for n = 1:numel(ucids.unit_cell_ids) - 1
        ism = ismember(G.cells.indexMap, ucids.unit_cell_ids{n});
        G.cells.unit(ism) = n;
        if n <= max(ucids.id.units) % units
            G.cells.Vcl(ism)     = ucids.Vcl(n);
            G.cells.juxt(ism, :) = repmat([0, 0], sum(ism), 1);
        else % fault cells
            id                   = n - max(ucids.id.units);
            G.cells.juxt(ism, :) = repmat(ucids.fault.juxt(id, :), ...
                                          sum(ism), 1); 
        end
    end
else
    for n = ucids.id.units
        ism = ismember(G.cells.indexMap, ucids.unit_cell_ids{n});
        G.cells.unit(ism) = n;
        G.cells.Vcl(ism)  = ucids.Vcl(n);
    end
    ism = ismember(G.cells.indexMap, ucids.unit_cell_ids{end});
    G.cells.unit(ism) = ucids.id.fault;
end
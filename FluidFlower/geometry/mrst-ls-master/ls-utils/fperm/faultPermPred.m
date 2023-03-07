function [kPred, throw, ContactsNodes, fcells, kFault0] = faultPermPred(G, ucids, name, optSmooth, optContacts, tlim)
%
%
%

%
if nargin < 6
    tlim.bot     = [];
    tlim.top.val = 0;
    tlim.top.z   = 0;
elseif ~isfield(tlim, 'bot')
    tlim.bot = [];
elseif ~isfield(tlim, 'top')
    tlim.top.val = 0;
    tlim.top.z   = 0;
end

% Preliminaries
ncf = diff(G.cells.facePos(1:2));                                           % number of faces per cell
fid = ucids.id.fault;
lid = ucids.id.faulted_units;
fc  = ucids.unit_cell_ids{fid};                                             % fault cells (all grid)
fcells = reshape(fc, numel(fc), 1);

%% Throw
% 1. Get faces bounding each layer (top, bottom) at each side of the fault
if G.griddim == 3
    %nc = 1:G.layerSize;
    fc = fc(fc < G.layerSize);
    nfc = numel(fc); 
    switch G.type{1}
        case 'triangleGrid'
            fmm=zeros(numel(lid), 2);
            [cwf, cwf2] = deal([], []);
            for n = 1:numel(lid)
                % layer cells
                lc  = [ucids.unit_cell_ids{lid(n)}];
                lc  = lc(lc < G.layerSize);
                
                % common face
                toadd = 0:ncf-1;
                faceid1 = repmat(G.cells.facePos(fc), 1, ncf);
                faceid2 = repmat(G.cells.facePos(lc), 1, ncf);
                faceid1 = bsxfun(@plus, faceid1, toadd); faceid2 = bsxfun(@plus, faceid2, toadd);
                cell_faces1 = G.cells.faces(faceid1); cell_faces2 = G.cells.faces(faceid2);
                comFace = intersect(cell_faces1, cell_faces2);
                
                [~, cmax] = max(G.faces.centroids(comFace, G.griddim));     % Hangingwall
                [~, cmin] = min(G.faces.centroids(comFace, G.griddim));     % Footwall
                fmm(n, :) = [comFace(cmax), comFace(cmin)];
                
                % cells (not fault) with common faces (all)
                nbors = G.faces.neighbors(comFace, :);
                cwf2   = [cwf2; nbors(~ismember(nbors, fc))];
                
                % cells (not fault) in HW -- for faults aligned with x dim
                nbors = G.faces.neighbors(comFace, :);
                [~, cid] = max([G.cells.centroids(nbors(:, 1), 2), ...
                                G.cells.centroids(nbors(:, 2), 2)], [], 2);
                cid   = diag(nbors(:, cid));
                cwf   = [cwf; cid(~ismember(cid, fc))];
                
            end
        otherwise                                                           % Hexes - cube / rectangular cuboid
            
    end
    
else                                                                        
    
end

% Check - plot
% uc = [ucids.unit_cell_ids{9:11}]; uc = uc(uc < G.layerSize);
% plotCellData(G, rock.poro, uc, 'faceAlpha', 0)
% hold on
% plotFaces(G, fmm(:,1),'FaceColor','r','FaceAlpha',1);
% plotFaces(G, fmm(:,2),'FaceColor','m','FaceAlpha',1);
% axis off

% 2. Compute throw at top and bottom of each offset layer
n2fn = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)';
zmm = zeros(size(fmm));
nl = numel(fmm(:,1));
for n = 1:nl % get depths at top nodes (FW) and bottom nodes of faces
    nod  = G.faces.nodes(ismember(n2fn, fmm(n,1)));
    nod2  = G.faces.nodes(ismember(n2fn, fmm(n,2)));
    zmm(n, :) = [max(G.nodes.coords(nod, 3)), min(G.nodes.coords(nod2, 3))];
end
t_lb = zmm(2:nl, 1) - zmm(1:nl-1, 2);  % throws at layer boundaries
if isempty(tlim.bot)
    tlim.bot = max(t_lb);                 % throw at fault bottom
end
t_lb = [tlim.bot; t_lb; tlim.top.val];              % add throw at fault bottom and top
z_t  = [zmm(1,1); zmm(1:nl-1, 2); tlim.top.z];      % z coordinates corresponding to throw (FW), bot to top

% 3. Throw of each fault cell based on linear interpolation
tSect = interp1(z_t, t_lb, G.cells.centroids(fc, G.griddim), 'linear');     % 2D section
throw = interp1(z_t, t_lb, G.cells.centroids(fcells, G.griddim), 'linear'); % ND
                                         
                                         
% get faces of contacts in both HW and FW
if optContacts == 1
    nall = 1:G.nodes.num;
    ncont = numel(zmm(:,1)) - 1;
    ContactsNodes = cell(ncont, 1);
    for c = 1:numel(zmm(:,1))-1
        hwn = nall(G.nodes.coords(:,3) == zmm(c+1, 1))';
        fwn = nall(G.nodes.coords(:,3) == zmm(c, 2))';
        ContactsNodes{c} = [fwn, hwn];
    end 
end

% Check - plot
% colormap(jet)
% plotCellData(G, throw, ucids.unit_cell_ids{fid})
% c = colorbar;
% hold on
% plotFaces(G, fmm(:,1),'FaceColor','r','FaceAlpha',1);
% plotFaces(G, fmm(:,2),'FaceColor','m','FaceAlpha',1);
% axis off
% t_lb

%% k Predictor
% Preliminaries: layer thickness in the HW
%lthw    = [abs(diff(zmm(:, 1))); zmm(end, 1)];                             % bot to top

% 1. For each fault cell, get cell in HW that has z closest to z of each
%    fault cell + corresponding throw
zhw     = G.cells.centroids(fc, G.griddim) + tSect;
dist    = pdist2(zhw, G.cells.centroids(cwf, G.griddim));
[~, id] = min(dist, [], 2);
cb      = cwf(id);                         % bottom cell in HW

% 2. Get all cells (in HW) offset by each cell in the fault
idc = repmat(cwf, 1, nfc);                                                  % size = N cells with shared face x N fault cells
id1 = repmat(G.cells.centroids(cwf, 3), 1, nfc);
idc1 = (id1 > G.cells.centroids(fc, 3)');
idc2 = (id1 < G.cells.centroids(cb, 3)');
ids = (idc1 + idc2) > 1;                                                    % ids of cells with shared face offset at each fault cell
%idc = idc(ids);     ids = cumsum(sum(ids, 1));

% 3. For each fault cell, find which units the offset cells belong to.
%    Then, compute thicknesses of each displaced unit and get their clay
%    fractions.
frac_dz  = cell(nfc, 1);
Vcl      = cell(nfc, 1);
%kFault0  = zeros(nfc,numel(rock.perm(1,:)));
for j  = 1:nfc
   tcells = [cb(j); idc(ids(:, j), j)];                                     % bottom cell must be added for all
   [ucells, id1] = unique(G.cells.unit(tcells), 'stable');
   nuc = numel(ucells);
   if nuc == 1
%        if tSect(j) == 0                                                   % fault throw = 0, perm = neighboring cell in HW (= FW)
%            kFault0(j, :) = rock.perm(tcells, :);
%        end                                                               
       frac_dz{j} = 1;
       Vcl{j}     = ucids.Vcl(ucells);
   else
       id2 = id1(2:end)-1;
       id  = zeros(numel(id1)+numel(id2), 1);
       id(1:2:end) = id1;   id(2:2:end) = id2;
       id  = [id; numel(tcells)];
       assert(numel(id)/2 == nuc);
       id  = reshape(id, 2, nuc);
       for k = 1:nuc                                                        % loop over each unit with offset cells
          iid = id(1,k):id(2,k);
          frac_dz{j}(k) = max(G.cells.centroids(tcells(iid), 3)) - ...
                          min(G.cells.centroids(tcells(iid), 3));
          if frac_dz{j}(k) == 0                                             % only one cell of given unit is offset
              if k == 1
                  iid0 = id(1,k+1):id(2,k+1);
              else
                  iid0 = id(1,k-1):id(2,k-1);
              end
              [~, idmin] = min(pdist2(G.cells.centroids(tcells(iid0), :), ...
                                   G.cells.centroids(tcells(iid), :)));
              iidf = iid0(idmin);
              frac_dz{j}(k) = abs(G.cells.centroids(tcells(iidf), 3) - ...
                                  G.cells.centroids(tcells(iid), 3));
          end
          Vcl{j}   = ucids.Vcl(ucells);
       end
       frac_dz{j} = frac_dz{j}./sum(frac_dz{j});
   end
   %dz{j} = frac_dz{j}.*tSect(j);
end

% 4. Compute SGR (required for all cases, since we need to provide props
%                 also where smears are discontinuous)
products = cellfun(@times, Vcl, frac_dz, 'UniformOutput', false);
nuof     = cellfun(@numel,frac_dz);
nuof     = abs(nuof - max(nuof));
for n = 1:nfc
   products{n} = [products{n} zeros(1, nuof(n))]; 
end
SGRSect =  sum(cell2mat(products), 2);                                      % raw as computed from cells
if optSmooth == true
    SGRSect = smooth(G.cells.centroids(fc, 3) , SGRSect, 0.03, 'rloess');   % Smooth SGR
end
nd = 5;
dpts = {round(G.cells.centroids(fc, G.griddim), nd), ...
        round(G.cells.centroids(ucids.unit_cell_ids{fid}, G.griddim), nd)}; % Round depths (in m) to nd decimals to ensure no difference in same depth due to machine error.
SGR  = interp1(dpts{1}, SGRSect, dpts{2}, 'linear');                        % Full SGR (extruded)
% Plot
%vals = ~isnan(SGR);
%plotCellData(G, SGR(vals), ucids.unit_cell_ids{fid}(vals))
%hold on
%plotCellData(G, rock.poro, cwf2(G.cells.unit(cwf2) == 5), 'FaceAlpha', 0, 'EdgeColor', 'r')
%set(gca,'Xdir','reverse')
%c = colorbar;

% 5. Compute predictor
kPred.SGR = SGR;
switch name
    case 'SGR'
        % already computed
        % kPred.vals = SGR; % not needed here since SGR is always passed.
        kPred.name = name;
    case 'SSF'
        error('SSF not coded yet')
    otherwise
        error('Predictor not supported. Type "help getFaultKPred"')
end

end


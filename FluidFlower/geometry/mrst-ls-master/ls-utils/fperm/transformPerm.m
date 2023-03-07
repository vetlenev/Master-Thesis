function [rock, fault] = transformPerm(G, rock, fault, opt)
%
%
%

% 1. Compute fault dip
fn_yz           = sort(fault.fn.fnodcoord(:,2:3), 1);
dfn_yz          = diff(fn_yz);
dip_p           = atan2d(dfn_yz(:,2), dfn_yz(:,1));                         % dips based on node linear conn.
y_cent          = fn_yz(1:end-1,2) + diff(fn_yz(:,2))/2;
[xData, yData]  = prepareCurveData( y_cent, dip_p );
dip_fcn         = fit( xData, yData, 'smoothingspline');                    % spline is the method used in Trelis for the nodeset

% 2. Compute theta (to rotate to global axis); 
%    Doc in https://petrowiki.org/Diagonalizing_the_permeability_tensor
fc        = fault.fcells;
fault.dip = dip_fcn(G.cells.centroids(fc, 3));
theta     = 90 - fault.dip;
% figure(83); plot(dip_cells, G.cells.centroids(unit_cells{end}, 3), '+b'); 
% set(gca,'Ydir','reverse'); grid on

% 3. Rotate
% clockwise rotation about the x' axis, since the along (in hz plane) axis 
% of the fault is aligned with the global x axis.
if isfield(fault.k, 'vals')
    kvals = fault.k.vals;
else
    kvals = fault.k;
end
for n=1:numel(fc)
    switch opt
        case 'rotateAboutX'
            T=[1 0 0; 0 cosd(theta(n)) sind(theta(n)); ...
               0 -sind(theta(n)) cosd(theta(n))];
            k_glob = T'*diag([kvals(n, 2), kvals(n, 1), kvals(n, 3)])*T;
    end
    assert(all(eig(k_glob) > 0))                                            % must be positive definite
    rock.perm(fc(n), 1) = k_glob(1,1);
    rock.perm(fc(n), 2) = k_glob(1,2);
    rock.perm(fc(n), 3) = k_glob(1,3);
    rock.perm(fc(n), 4) = k_glob(2,2);
    rock.perm(fc(n), 5) = k_glob(2,3);
    rock.perm(fc(n), 6) = k_glob(3,3);
end

end
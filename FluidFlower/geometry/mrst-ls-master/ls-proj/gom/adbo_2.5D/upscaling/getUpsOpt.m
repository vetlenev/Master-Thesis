function opt = getUpsOpt(opt)
%
%
%
if strcmp(opt.fault, 'predict') || strcmp(opt.fault, 'predict_3D') 
    opt.sg = 'sandClay';
    if strcmp(opt.strati, 'GoM') 
        opt.maxPerm = 175; % mD, max perm of Amp B interval (sand layers)
        
        pf = [12652.3 12626.6 12600.9 12575.2 12549.5; ...
              1751.3 1726.3 1701.3 1676.3 1651.3]; % substitute for each throw window
        ph = 0;
        lx = pf(1,1:end-1) - pf(1,2:end);
        lz = pf(2,1:end-1) - pf(2,2:end);
        tap = sqrt(lx.^2 + lz.^2);
        zmax = pf(2,1:end-1) + diff(pf(2,:))/2;
        
        pdip = [12623.1 12685.6; 1700 1689.24]; % substitute for each throw window
        lx = pdip(1,2) - pdip(1,1);
        lz = pdip(2,1) - pdip(2,2);
        dip_lyr = atand(lz/lx);
        
        lx = pf(1,1) - pf(1,end);
        lz = pf(2,1) - pf(2,end);
        fdip = atand(lz/lx);
        
        if strcmp(opt.window, 'famp1') % bottom throw window
            opt.thick = {[115.6143 28.8949], [37.6113 37.6861 37.6113 31.6005]};       % from mesh
            if opt.clay_vcl == 0.7
                opt.vcl = {[0.3 0.7], [0.2878,0.7,0.2878,0.7]};
            elseif opt.clay_vcl == 0.65
                opt.vcl = {[0.3 0.65], [0.2878,0.65,0.2878,0.65]};
            elseif opt.clay_vcl == 0.55
                opt.vcl = {[0.3 0.55], [0.2878,0.55,0.2878,0.55]};
            elseif opt.clay_vcl == 0.45
                opt.vcl = {[0.3 0.45], [0.2878,0.45,0.2878,0.45]};
            end
            opt.dip = [0, -12.0136];     % from mesh
            opt.fDip = 41.6345;         % from mesh
            opt.zf = [200, 200];        % 200m
            opt.zmax = {[1912 1861], [1934 1909 1884 1860]};        % center of each layer, from mesh
        
        elseif strcmp(opt.window, 'famp2')
            opt.thick = {[36.9255 35.8537 36.8537 36.3111], [36.5042 36.5042 36.4314 36.5042]};  
            if opt.clay_vcl == 0.7
                opt.vcl = {[0.2878, 0.7, 0.2878, 0.7], [0.2878, 0.7, 0.2878, 0.7]};
            elseif opt.clay_vcl == 0.65
                opt.vcl = {[0.2878, 0.65, 0.2878, 0.65], [0.2878, 0.65, 0.2878, 0.65]}; 
            end
            opt.dip = [0, -13.8951];     
            opt.fDip = 43.2508;        
            opt.zf = [200, 200];        
            opt.zmax = {[1837.5 1812.5 1787.5 1762.5], [1837.5 1812.5 1787.5 1762.5]};
            
        elseif strcmp(opt.window, 'famp3')
            opt.thick = {[35.8537, 35.8537, 35.8537, 35.8537], [35.8537, 35.8537, 35.8537, 35.8537]};  
            if opt.clay_vcl == 0.7
                opt.vcl = {[0.2878, 0.7, 0.2878, 0.7], [0.2878, 0.7, 0.2878, 0.7]}; 
            end
            opt.dip = [0, -9.7683];     
            opt.fDip = 43.8;        
            opt.zf = [200, 200];        
            opt.zmax = {[1738.8 1713.8 1688.8 1663.8], [1738.8 1713.8 1688.8 1663.8]};
            
        elseif strcmp(opt.window, 'famp4')
            opt.thick = {[35.8537, 35.8537, 35.8537, 35.9255], [35.8537, 35.8537, 35.8537, 35.9255]};  
            if opt.clay_vcl == 0.7
                opt.vcl = {[0.2878, 0.7, 0.2878, 0.7], [0.2878, 0.7, 0.2878, 0.7]}; 
            end
            opt.dip = [0, -4.9456];     
            opt.fDip = 44.1811;        
            opt.zf = [200, 200];        
            opt.zmax = {[1638.8 1613.8 1588.8 1563.8], [1638.8 1613.8 1588.8 1563.8]};
            
        elseif strcmp(opt.window, 'famp5')
            opt.thick = {[35.8537 35.8537 35.8537 35.8537], [37.4901 35.2847 35.3553 35.2847]};  
            if opt.clay_vcl == 0.7
                opt.vcl = {[0.2878, 0.7, 0.2878, 0.7], [0.2878, 0.7, 0.2878, 0.7]};
            end
            opt.dip = [0, -5.2221];     
            opt.fDip = 45.0685;        
            opt.zf = [200, 200];        
            opt.zmax = {[1538.8 1513.82 1488.75 1463.99], [1538.8 1513.82 1488.75 1463.99]};
            
        elseif strcmp(opt.window, 'famp6')
            opt.thick = {[28.2932, 33.1042, 33.1699, 33.1042], 127.6715};  
            if opt.clay_vcl == 0.7
                opt.vcl = {[0.2878, 0.7, 0.2878, 0.7], 0.2}; 
            end
            opt.dip = [0, -5];     
            opt.fDip = 46.0685;        
            opt.zf = [200, 200];        
            opt.zmax = {[1440.6 1417.5 1392.5 1367.5], 1400};
            
        elseif strcmp(opt.window, 'smearCross')
            opt.thick = {[85 15], [100]};       
            opt.vcl = {[0.2 0.8], [0.3]};
            opt.dip = [0, 0];     
            opt.fDip = 60;          
            opt.zf = [200, 200];         
            opt.zmax = {[1912 1861], [1900]};
            
        end
        
    end
    
% x = [12652.3 12678];
% z = [1751.3 1776.3];
% beta = atand(diff(z)/diff(x));
% tap = diff(x)/cosd(beta)
% 
% f=figure(1);
% plotToolbar(G, log10(rock.perm(:,1)/(milli*darcy)), 1:G.layerSize);
% view([-90, 0])
% set(gca,'Ydir','reverse')
% axis equal
% ids = ismember(ucids.unit_cell_ids{27}, 1:G.layerSize);
% plotGrid(G, ucids.unit_cell_ids{27}(ids), 'faceColor', 'none', 'edgecolor', 'r');
% ylim([11000 16000])
% zlim([1300 2100])
% colormap(copper)
% f.Position = [100 100 1500 400];

elseif strcmp(opt.fault, 'test')
    opt.sg = 'sandClay';
end

end
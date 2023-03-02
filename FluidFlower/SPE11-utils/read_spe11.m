filename = 'draft_spe11_with_facies_markers.geo';

[pts, loops, facies] = parse_spe11_geo(filename);


% plot each faciles in a different color
colors = {[0 0.4470 0.7410], 
          [0.8500 0.3250 0.0980],
          [0.9290 0.6940 0.1250],
          [0.4940 0.1840 0.5560],
          [0.4660 0.6740 0.1880],
          [0.3010 0.7450 0.9330],
          [0.6350 0.0780 0.1840]};

% plot result
figure;
for f_ix = 1:numel(facies)
    facie = facies{f_ix};
    
    for poly_ix = facie'
        
        poly = pts(loops{poly_ix}, :);
        poly = [poly; poly(1,:)];
        patch(poly(:,1), poly(:,2), colors{f_ix}); hold on
    end
end

for f_ix = 1:numel(facies)
    facie = facies{f_ix};   
    for poly_ix = facie'        
        poly = pts(loops{poly_ix}, :);
        poly = [poly; poly(1,:)];       
        text(mean(poly(:,1)), mean(poly(:,2)), string(poly_ix)); hold on
    end
end
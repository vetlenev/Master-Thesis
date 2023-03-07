%
%
%

% changes vs each model
%          k pc kr Im D
nsim = [20, 10, 10];
changes = [11 4 1 4 5;
           7 8 0 4 1;
           10 6 0 6 1];
          

sg_area = [91 3.22 89];                % [inf t=54, med t=154, sup t=154]  cm2
giw_area = [234.7 33.21 251.53];
sg_area_models = [80.1 2.71 86.5;
                  81.24 2.71 88.7;
                  90.64 2.71 102.61]/1.05;
giw_area_models = [263.2 34.7 228.9;
                   255.49 38.04 272.62;
                   244.48 34.963 259.32]/1.05;

t_finger = [315 420];                   % min, from experiment
t_finger_models = [365 395;
                   365 395;
                   355 405];            % "

close all
clrs = hot(10);
labls = {'Model 1', 'Model 2', 'Model 3'};
latx = {'Interpreter', 'latex'};
tit={'Gas in Finf', 'Brine with gas in Finf', 'Gas in Fmid', ...
     'Brine with gas in Fmid','Gas in Fsup', 'Brine with gas in Fsup'};
for n=1:3
   f1 = figure(1);
   hold on
   b1 = bar(n, nsim(n), 'k', 'BarWidth', 0.5);
   hold off
   grid on
   xlim([0.5 3.5])
   xticks([1 2 3])
   title('Number of simulations to match (4mm mesh)')
    
   f3 = figure(3);
   subplot(1,6,1)
   hold on
   b = bar(n, abs(sg_area_models(n, 1)-sg_area(1)));
   b.FaceColor = 'flat';
   b.CData = clrs(n,:);
   grid on
   hold off
   subplot(1,6,2)
   hold on
   b = bar(n, abs(giw_area_models(n, 1)-giw_area(1)));
   b.FaceColor = 'flat';
   b.CData = clrs(n,:);
   grid on
   hold off
   subplot(1,6,3)
   hold on
   b = bar(n, abs(sg_area_models(n, 2)-sg_area(2)));
   b.FaceColor = 'flat';
   b.CData = clrs(n,:);
   grid on
   hold off
   subplot(1,6,4)
   hold on
   b = bar(n, abs(giw_area_models(n, 2)-giw_area(2)));
   b.FaceColor = 'flat';
   b.CData = clrs(n,:);
   subplot(1,6,5)
   hold on
   b = bar(n, abs(sg_area_models(n, 3)-sg_area(3)));
   b.FaceColor = 'flat';
   b.CData = clrs(n,:);
   subplot(1,6,6)
   hold on
   b = bar(n, abs(giw_area_models(n, 3)-giw_area(3)));
   b.FaceColor = 'flat';
   b.CData = clrs(n,:);
   if n==3
        subplot(1,6,1)
        xlabel('Simulation nr.', 'fontsize', 12)
        ylabel('$\vert A_\mathrm{model} - A \vert$ [cm$^2$]', latx{:}, 'fontsize', 12)
        grid on
        title(tit{1}, 'fontSize', 12)
        subplot(1,6,2)
        grid on
        title(tit{2}, 'fontSize', 12)
        subplot(1,6,3)
        grid on
        title(tit{3}, 'fontSize', 12)
        subplot(1,6,4)
        grid on
        title(tit{4}, 'fontSize', 12)
        subplot(1,6,5)
        grid on
        title(tit{5}, 'fontSize', 12)
        subplot(1,6,6)
        grid on
        title(tit{6}, 'fontSize', 12)
   end
    
   f4 = figure(4);
   subplot(1,2,1)
   hold on
   b = bar(n, abs(t_finger_models(n, 1)-t_finger(1)));
   b.FaceColor = 'flat';
   b.CData = clrs(n,:);
   subplot(1,2,2)
   hold on
   b = bar(n, abs(t_finger_models(n, 2)-t_finger(2)));
   b.FaceColor = 'flat';
   b.CData = clrs(n,:);
   if n == 3
       subplot(1,2,1)
       grid on
       title('Finger Finf', 'fontSize', 11)
       ylabel('$\vert t_\mathrm{model} - t \vert$ [min]', latx{:}, 'fontsize', 12)
       xlabel('Model')
       xticks([1 2 3])
       
       subplot(1,2,2)
       grid on
       title('Finger Fsup', 'fontSize', 11)
       xticks([1 2 3])
   end
   
   f2 = figure(2);
   hold on
   b = bar(n, changes(n,:), 'stacked');
   ax = gca;
   b(1).FaceColor = clrs(2,:);
   b(2).FaceColor = clrs(3,:);
   b(3).FaceColor = clrs(5,:);
   b(4).FaceColor = clrs(7,:);
   b(5).FaceColor = clrs(9,:);
   hold off
   grid on
   if n==3
      xt = get(ax, 'XTick');
      set(gca, 'XTick', xt([3 5 7]), 'XTickLabel', labls)
      ylabel('Nr of changes', 'fontSize', 12)
   end
end
yd = get(b, 'YData');
ylb = {'$k$' '$P_c$' '$k_r$' '$I$' '$D$'};
barbase = cumsum([zeros(size(changes,1),1) changes(:,1:end-1)],2);
ylblpos = changes/2 + barbase;
text(xt(3)*ones(1,size(changes,2)), ylblpos(1,:), ylb, ...
    'HorizontalAlignment','center', 'color', 'k', latx{:}, ...
    'fontSize', 12)
f2.Position = [50, 50, 400, 300];

f1.Position = [50 50 325, 300];

f3.Position = [100, 100, 1200, 250];

f4.Position = [50 50 400, 300];
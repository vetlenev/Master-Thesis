% This is for model 1!
sg_area = [91 3.22 89];                % [inf t=54, med t=154, sup t=154]  cm2
giw_area = [234.7 33.21 251.53];       % ["]                               cm2

sg_v_model = [16.59 2.56 14.3;             % % first to last 4mm folderName based on date modified in GRS3
                 16.41 2.56 14.16; 
                 1.70 0 0;
                 123.49 23.16 87.77;
                 100.97 nan nan;
                 91.64 8.16 31.81;
                 80.99 33.28 39.42;
                 82.57 2.71 90.03;
                 88.47 2.71 91.19;
                 80.86 2.71 90.34;
                 80.27 2.71 89.14;
                 80.11 2.82 88.84;
                 80.11 2.82 69.89;
                 80.11 2.71 87.85;
                 80.73 2.71 88.11;
                 80.73 2.71 86.62;
                 80.27 2.71 86.62;
                 80.11 2.71 86.47;
                 80.11 2.71 86.69;
                 80.1 2.71 86.5];  
             
giw_v_model = [322.44 66.22 315.42;
                  323.33 66.27 316.36;
                  345 415.55 126.63;
                  312.24 87.43 275.37;
                  286.49 nan nan;
                  251.59 19.6 85.36;
                  253.02 105.19 159.1;
                  255.74 24.28 229.7;
                  252.05 23.78 229.47;
                  254.58 24.65 233.03;
                  260.28 25.16 230.59;
                  263.29 26.45 231.42;
                  264 29.75 231.26;
                  263.54 28.13 229.63;
                  259.62 27.44 230.53;
                  259.78 31.76 230.47;
                  260.61 31.43 230.81;
                  264.21 32.6 230.5;
                  283.82 31.47 226.19;
                  263.2 34.7 228.9];   % "

sg_area_model = sg_v_model./1.05;
giw_area_model = giw_v_model./1.05;
              
nsim = size(sg_area_model, 1);
latx = {'Interpreter', 'latex'};
fontsz = {'fontSize', 14};
tit={'Gas in Finf', 'Brine with gas in Finf', 'Gas in Fmid', ...
     'Brine with gas in Fmid','Gas in Fsup', 'Brine with gas in Fsup'};
f = figure(11);
colrs = colormap(hot(nsim));
for n=1:nsim
    subplot(1,6,1)
    hold on
    b = bar(n, abs(sg_area_model(n, 1)-sg_area(1)));
    b.FaceColor = 'flat';
    b.CData = colrs(n,:);
    subplot(1,6,2)
    hold on
    b = bar(n, abs(giw_area_model(n, 1)-giw_area(1)));
    b.FaceColor = 'flat';
    b.CData = colrs(n,:);
    subplot(1,6,3)
    hold on
    b = bar(n, abs(sg_area_model(n, 2)-sg_area(2)));
    b.FaceColor = 'flat';
    b.CData = colrs(n,:);
    subplot(1,6,4)
    hold on
    b = bar(n, abs(giw_area_model(n, 2)-giw_area(2)));
    b.FaceColor = 'flat';
    b.CData = colrs(n,:);
    subplot(1,6,5)
    hold on
    b = bar(n, abs(sg_area_model(n, 3)-sg_area(3)));
    b.FaceColor = 'flat';
    b.CData = colrs(n,:);
    subplot(1,6,6)
    hold on
    b = bar(n, abs(giw_area_model(n, 3)-giw_area(3)));
    b.FaceColor = 'flat';
    b.CData = colrs(n,:);
    if n==nsim
        subplot(1,6,1)
        xlabel('Simulation nr.', fontsz{:})
        ylabel('$\vert A_\mathrm{model} - A \vert$ [cm$^2$]', latx{:}, fontsz{:})
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
end
f.Position = [0, 0, 1500, 250];
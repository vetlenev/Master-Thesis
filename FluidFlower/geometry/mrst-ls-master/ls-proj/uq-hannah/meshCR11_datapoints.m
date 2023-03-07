%% Make datapoints for mesh conforming to boundary, fault, folded layers

% Boundary
stratiPoints = struct('boundary', [], 'lines', []);
stratiPoints.boundary = [0, 500; ...
                         2000, 500; ...
                         2000, 2500; ...
                         0, 2500; ...
                         0, 500];
                     
fT = 10;     % fault thickness (m)

% Layer lines (top to bottom)
% Hangingwall
stratiPoints.lines{3} = [0, 1700; ...
                         464.7346, 1700]; 
stratiPoints.lines{4} = [0, 1550; ...
                         491.1837, 1550];
stratiPoints.lines{5} = [0, 1450; ...
                         508.8163, 1450];
stratiPoints.lines{6} = [0, 1300; ...
                         535.2654, 1300];
% Footwall             
stratiPoints.lines{7} = [2000, 1575; ...
                         500+fT-tand(10)*75, 1575]; 
stratiPoints.lines{8} = [2000, 1425; ...
                         500+fT+tand(10)*75, 1425];
stratiPoints.lines{9} = [2000, 1325; ...
                         500+fT+tand(10)*175, 1325];
stratiPoints.lines{10} = [2000, 1175; ...
                          500+fT+tand(10)*325, 1175];

% Fault lines (There is 1 fault-like structure, but we need two lines to 
%              define its thickness)
stratiPoints.lines{1} = [323.6730, 2500;
                         stratiPoints.lines{3}(2,1), 1700; 
                         500-tand(10)*175, 1675;
                         stratiPoints.lines{4}(2,1), 1550; 
                         stratiPoints.lines{5}(2,1), 1450; 
                         500+tand(10)*75, 1425;
                         stratiPoints.lines{6}(2,1), 1300;
                         500+tand(10)*325, 1175;
                         676.3270, 500];
stratiPoints.lines{2} = [stratiPoints.lines{1}(1)+fT, 2500;
                         500+fT-tand(10)*175, 1675;
                         stratiPoints.lines{7}(2,1), 1575;
                         500+fT-tand(10)*50, 1550;
                         stratiPoints.lines{8}(2,1), 1425;                         
                         stratiPoints.lines{9}(2,1), 1325;
                         500+fT+tand(10)*200, 1300;
                         stratiPoints.lines{10}(2,1), 1175;
                         stratiPoints.lines{1}(end,1)+fT, 500];   
                     
stratiPoints.lines{11} = [stratiPoints.lines{1}(3,:); stratiPoints.lines{2}(2,:)];
stratiPoints.lines{12} = [stratiPoints.lines{1}(4,:); stratiPoints.lines{2}(4,:)];
stratiPoints.lines{13} = [stratiPoints.lines{1}(6,:); stratiPoints.lines{2}(5,:)];
stratiPoints.lines{14} = [stratiPoints.lines{1}(7,:); stratiPoints.lines{2}(7,:)];
stratiPoints.lines{15} = [stratiPoints.lines{1}(8,:); stratiPoints.lines{2}(8,:)];
                    
% Save data
% Run the following command to see directory where datapoints will be saved
% >> pwd 
save('meshCR11_datapoints.mat', 'stratiPoints');
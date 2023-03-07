%% Make datapoints for mesh conforming to boundary, fault, folded layers
%
% We draw the stratigraphy using paths in inkscape. Then, we export nodes
% by selecting a given path and then doing extensions>export>export XY. Finally,
% we transform the y values using the actual y location (in mm) of the 
% leftmost node, since they are initially not absolute coordinates that we
% want. The extension export XY is external and was downloaded from here
% https://alpha.inkscape.org/vectors/www.inkscapeforum.com/viewtopicdb3f.html?t=8826
% and added to inkscape.
%

% Boundary
h = 47; % height (cm)
l = 89.7; % length
to_m = 1e-2; % convert measuring/drawing units (cm) to m
stratiPoints = struct('boundary', [], 'lines', [], 'wells', []);
stratiPoints.boundary = to_m*[0, 0; ...
                              l, 0; ...
                              l, -h; ...
                              0, -h; ...
                              0, 0];

% transforms
yn = @(y) -h + y;
%                                       %top %actual y coord leftmost node
to_mesh = @(vals, yl) to_m*[vals(:,1) -cumsum([h-yl; diff(vals(:,2))])];
                          
% Wells (to have cells with correct centroids)
stratiPoints.wells{1} = to_m*[45, yn(2)];
stratiPoints.wells{2} = to_m*[16.5, yn(16)];
stratiPoints.wells{3} = to_m*[8, yn(8)];
stratiPoints.wells{4} = to_m*[81.7, yn(8)];
stratiPoints.wells{5} = to_m*[7.9, yn(16)];
stratiPoints.wells{6} = to_m*[24.3, yn(24.3)];
stratiPoints.wells{7} = to_m*[2, yn(2)];
stratiPoints.wells{8} = to_m*[87.7, yn(2)];

% Fault lines (1 fault, 2 lines are shear zone boundaries)
stratiPoints.lines{1} = to_mesh([9.638,296.99999
                                 10.181953,293.00334
                                 10.897483,287.76403
                                 11.095287,286.22501
                                 11.879124,280.4356
                                 12.339991,276.8547
                                 12.567708,275.26636], 0);
                          
stratiPoints.lines{2} = to_mesh([10.317505,296.99999
                                 11.097145,291.2205
                                 11.817121,285.80632
                                 12.024841,284.2855
                                 12.849296,278.32276
                                 13.247213,275.26636], 0);

stratiPoints.lines{3} = [stratiPoints.lines{1}(end,1), stratiPoints.lines{1}(end,2); ...
                         stratiPoints.lines{2}(end,1), stratiPoints.lines{2}(end,2)]; % fault closure

% Layer lines (top to bottom and left to right)
% HW
stratiPoints.lines{4} = to_mesh([10.181953,293.00334
                                8.7493535,292.61374
                                6.7435675,292.42271
                                4.8127065,292.44631
                                2.8924692,292.58503
                                0.74284868,292.81341
                                0.0,293.0336], stratiPoints.lines{1}(2,2)/to_m + h);
           
stratiPoints.lines{5} = to_mesh([10.897483,287.76403
                                        9.3822018,287.33899
                                        7.3764158,287.14796
                                        5.4455547,287.17156
                                        3.5253174,287.31028
                                        1.3756968,287.53866
                                        0.0,287.75899], stratiPoints.lines{1}(3,2)/to_m + h);

stratiPoints.lines{6} = to_mesh([11.095287,286.22501
                                        9.6006779,285.81179
                                        7.5948919,285.62076
                                        5.664031,285.64436
                                        3.7437936,285.78308
                                        1.5941731,286.01146
                                        0.0,286.23165], stratiPoints.lines{1}(4,2)/to_m + h);
                                    
stratiPoints.lines{7} = to_mesh([11.879124,280.4356
                                        10.229196,280.09897
                                        8.2234103,279.92465
                                        6.626636,279.94825
                                        5.1407115,280.07026
                                        2.9576825,280.26524
                                        1.3042999,280.44898
                                        0.0,280.71299], stratiPoints.lines{1}(5,2)/to_m + h);
                                    
stratiPoints.lines{8} = to_mesh([12.567708,275.26636
                                        6.1806079,275.26636
                                        3.9422255,275.35116
                                        1.9711127,275.5182
                                        0.33408691,275.75206
                                        0.0,275.85196], stratiPoints.lines{1}(end,2)/to_m + h);
% FW
stratiPoints.lines{9} = to_mesh([11.097145,291.2205
                                 14.566189,291.93129
                                 17.128874,292.19451
                                 20.398872,292.31599
                                 24.250656,292.22559
                                 28.753911,291.87599
                                 33.007786,291.25912
                                 38.620446,290.59094
                                 45.970359,290.19004
                                 52.919365,289.85595
                                 61.799107,289.6673
                                 73.327381,290.13977
                                 79.563988,290.61224
                                 85.517113,291.46269
                                 89.7,292.40763], ...
                                 stratiPoints.lines{2}(2,2)/to_m + h);

stratiPoints.lines{10} = to_mesh([11.817121,285.80599
                                  14.726121,286.40663
                                  17.288806,286.66985
                                  20.558804,286.79133
                                  24.410588,286.70093
                                  28.913843,286.35133
                                  33.167718,285.73446
                                  38.780379,285.06628
                                  46.130292,284.66538
                                  53.079298,284.33129
                                  61.959041,284.14264
                                  73.487315,284.61511
                                  79.72392,285.08758
                                  85.677045,285.93803
                                  89.7,286.88297], ...
                                  stratiPoints.lines{2}(3,2)/to_m + h);
                     
stratiPoints.lines{11} = to_mesh([12.024841,284.2855
                                  14.926922,284.87818
                                  17.584101,285.17683
                                  20.877722,285.2865
                                  24.811493,285.14683
                                  30.083148,284.66359
                                  34.337023,283.97991
                                  39.436567,283.34717
                                  46.330745,282.97764
                                  54.014742,282.54333
                                  63.575516,282.27975
                                  72.952779,282.79374
                                  80.05801,283.36644
                                  85.075692,284.21689
                                  89.7,285.29546], ...
                                  stratiPoints.lines{2}(4,2)/to_m + h);

stratiPoints.lines{12} = to_mesh([12.849296,278.32276
                                  15.705476,278.57118
                                  18.52802,278.68084
                                  21.721415,278.62344
                                  26.189725,278.51718
                                  31.093884,278.20098
                                  35.280942,277.48389
                                  40.380487,276.85115
                                  47.274663,276.48162
                                  54.958661,276.04731
                                  65.905896,275.76703
                                  73.8967,276.29772
                                  81.001929,276.87042
                                  86.019612,277.72087
                                  89.7,278.79944], ...
                                  stratiPoints.lines{2}(5,2)/to_m + h);

stratiPoints.lines{13} = to_mesh([13.247213,275.26636
                                  16.028223,275.5106
                                  18.850767,275.62026
                                  22.044162,275.56286
                                  26.512473,275.45659
                                  31.416632,275.14039
                                  35.60369,274.42331
                                  40.703234,273.79056
                                  47.597411,273.42104
                                  55.281408,272.98672
                                  66.228644,272.70644
                                  74.219447,273.23713
                                  81.324676,273.80983
                                  86.342359,274.66028
                                  89.7,275.73886], ...
                                  stratiPoints.lines{2}(6,2)/to_m + h);
                      
% Above
stratiPoints.lines{14} = to_mesh([0.0,267.06399
                                  4.9136905,267.44224
                                  10.110863,267.63123
                                  15.497024,267.91471
                                  22.206101,268.1982
                                  26.836309,268.2927
                                  32.411458,268.00922
                                  40.821427,267.06428
                                  50.270833,265.93035
                                  57.357885,265.36338
                                  64.255951,265.36338
                                  71.531992,265.83585
                                  82.87128,267.15877
                                  89.7,268.76517], 29.936);

stratiPoints.lines{15} = to_mesh([0.0,       262.0379
                                  5.0988895, 262.41583
                                  10.296062, 262.60482
                                  15.682224, 262.8883
                                  22.391301, 263.17179
                                  27.021508, 263.26629
                                  32.596657, 262.98281
                                  41.006626, 262.03787
                                  50.456032, 260.90394
                                  57.543084, 260.33697
                                  64.44115,  260.33697
                                  71.717192, 260.80944
                                  83.056479, 262.13236
                                  89.7,      263.73876], 34.962);


% Save data
save('mediumMesh_datapoints.mat', 'stratiPoints');
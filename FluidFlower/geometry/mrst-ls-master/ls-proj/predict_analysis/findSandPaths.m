%% Find if any object in a binary img extends L to R or Top to Bot
%
%

% Test matrix
dim = 5;
A = zeros(dim);
id1 = [3 4 9 10 15 11 12 17 18 24];
A(id1) = 1;

% Connectivity
conn2D = 8;     % MPFA, 4 for TPFA
L = labelmatrix(bwconncomp(A, conn2D));
B = regionprops(L, 'BoundingBox');
imshow(label2rgb(L,'jet'), 'InitialMagnification','fit');

% Find extension of each object
ext.x = arrayfun(@(s) s.BoundingBox(3), B)';
ext.z = arrayfun(@(s) s.BoundingBox(4), B)';
if any(ext.x == dim), disp('Connected in x dim (Left to Right)'), end
if any(ext.z == dim), disp('Connected in z dim (Top to Bot)'), end
function rand_pt = randomPointOnLine(p1, p2, rng_type)
% Create a random point on line segment from p1 to p2
%   Detailed explanation goes here
r = max(0, min(1, 0.2*rand(1) + 0.5));
if strcmp(rng_type, 'uniform') 
    rand_pt = r*p1 + (1-r)*p2;
elseif strcmp(rng_type, 'nonuniform')    
    rand_pt = sqrt(r)*p1 + (1-sqrt(r))*p2;
end

end


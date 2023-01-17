function [xp,yp,I,seglabel] = border_finderI(I,nt)
% This function finds the uniformly distributed coordinates of the 
% periphery of the embryo 
% 
% Input:
%       I: Image in uint16
%       nt: Final number of points
% Output: 
%       xp: x coordinates of the border
%       yp: y coordinates of the border
%       I: segmented image

[I, seglabel]   = segment_embryo_image(I);
Borders         = bwboundaries(seglabel,'noholes');
points          = Borders{1};
xp              = points(:,2); 
yp              = points(:,1);
N               = length(xp);
xp              = smooth_periodic(xp,round(0.05*N));       % made this 0.05% to reduce over-smoothing.
yp              = smooth_periodic(yp,round(0.05*N));
[xp,yp]         = makePointsUniform(xp,yp,nt);
[xp,yp]         = makePointsUniform(xp,yp,nt);


% check by plotting
%{
I = max(I,[],3); 
imshow(I,[]) 
hold on
plot(xp,yp,'*')
%}
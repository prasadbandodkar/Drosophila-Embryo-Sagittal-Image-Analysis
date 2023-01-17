function [xp_in,yp_in,xp_in2,yp_in2] = findNormals(xp,yp,varargin)
%Finds points locally normal to xp,yp, that are Yhatmax closer to xc,yc
%
%function [x_in,y_in] = find_normals(xp,yp,xc,yc,varargin)
%
% This function takes a set of points, (xp,yp), that are presumably the
% boundary of an embryo, and finds another set of points, (x_in,y_in) that
% are locally normal to the boundary, but are moved "in" by Yhatmax pixels.
%
% 
% Optional argument varargin can consist of the following things:
%	(1) "Yhatmax": the depth into the embryo that we keep. Default, 60 pxl.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	(2) "dthetamax": the maximum angle (in degrees) between two line
%		segments we allow before determining that the inner points need to
%		be re-written. Default, 45 degrees.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%


%load('Sample/sample_findnormals.mat')
%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	Yhatmax = varargin{iArg}; else
	Yhatmax = 60;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	dthetamax = varargin{iArg}; else
	dthetamax = 45;
end%, iArg = iArg + 1;

ns = length(xp)-1;
xp_in  = zeros(ns+1,1);   yp_in  = xp_in;
xp_in2 = xp_in;           yp_in2 = xp_in;


% First determine the vector that is tangent to the points.
a1  = xp(2:end) - [xp(end-1);xp(1:end-2)];
a2  = yp(2:end) - [yp(end-1);yp(1:end-2)];



%
% Then, determine the local normals. The inner tanget vector is (a1,a2,0)
% and the local normal is (b1,b2,0). The formulae are calculated such that
% the dot product of the two vectors is equal to zero and the cross product
% is equal to 1.
%
b2  = a1./(a1.^2 + a2.^2);
b1  = -(a2./a1).*b2;
d   = sqrt(b1.^2+b2.^2);
b1  = b1./d;
b2  = b2./d;

%xp_in(1:end) = xp(1:end) + Yhatmax*b1;
%yp_in(1:end) = yp(1:end) + Yhatmax*b2;

xp_in(1:end-1) = xp(1:end-1) + Yhatmax*b1;
yp_in(1:end-1) = yp(1:end-1) + Yhatmax*b2;
xp_in(end) = xp_in(1);
yp_in(end) = yp_in(1);

xp_in2(2:end) = xp(2:end) + Yhatmax*b1;
yp_in2(2:end) = yp(2:end) + Yhatmax*b2;
xp_in2(1) = xp_in(end);
yp_in2(1) = yp_in(end);


end



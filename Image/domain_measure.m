function [t,raw,s] = domain_measure(I,xp,yp,varargin)
%Finds domains of gene expression in a cross-sectioned embryo.
%
%function [t,raw,s] = domain_measure(I,xp,yp,varargin)
%
% This function takes a truecolor cross-sectional image of an embryo and
% finds the fluorescence intensity (gene expression) as a function of
% fractions of total circumference.
%
% Inputs:
%
% "I": can be either a uint8/uint16 grayscale image, or string of filename
% "xp,yp": the points along the perimeter of the embryo
% 
% Optional argument varargin can consist of the following things:
%	* "depthInsideImage":  the distance, in pixels into the embryo for
%       analysis
%		Default, 30 pixels.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "ns": choice for number of bins in "s" (pseudoarclength) when
%		measuring the intensity around the periphery.  Must be an even
%		number.  Default, 300.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% Outputs:
%
% "t": npts-by-nChannels array of smoothened fluorescent intensities.
% "raw": same as "t", but the raw (non-smoothened) data.
% "s": arclength, going from -1 to +1, with "npts" points.
% 
% PB Edit - 01/11/2020
% This function now finds the inner border, makes the border uniform (just 
% like the external border) and then finds the domains.

% load('Sample/domainMeasure.mat')

%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	depthInsideImage = varargin{iArg}; else
	depthInsideImage = 30;
end, iArg = iArg +1;
if nArg >= iArg && ~isempty(varargin{iArg})
	npts = varargin{iArg}; else
	npts = 350;
end

%
% Reading in "I"
%
if ischar(I)
	I = imread(I);
end
[m,n,o] = size(I);


% perimeter calculation
dx = diff(xp);              dy = diff(yp);  
ds = sqrt(dx.^2 + dy.^2);   perim = sum(ds);   
ss = 2*[0;cumsum(ds)]/perim - 1;
ss(end) = 1;

% find inner boundary
[xp_in,yp_in]  = findNormals(xp,yp,depthInsideImage);   
%{
[~,~,xc,yc,] = ellipse_fit(xp,yp);
[x_in,y_in] = find_normals(xp,yp,xc,yc);
%}


% check by plotting
%{
figure
imshow(I(:,:,2),[])
hold on
plot(xp,yp,'*')
hold on
plot(xp_in,yp_in,'*')
%plot(x_in,y_in,'*')
for i=1:length(xp)
    plot([xp(i),xp_in(i)],[yp(i),yp_in(i)])
    %plot([xp(i),x_in(i)],[yp(i),y_in(i)])
end
%}


%
% Building raw data
%
% defining image coordinates
x_im = (1:n);              y_im = (1:m)';
X    = repmat(x_im,m,1);   Y    = repmat(y_im,1,n);


% calculating variables for use in loop
deltax          = diff(xp);
deltay          = diff(yp);
deltax_across   = xp_in(1:end-1) - xp(1:end-1);
deltay_across   = yp_in(1:end-1) - yp(1:end-1);
deltax_across1  = xp_in(2:end) - xp(1:end-1);
deltay_across1  = yp_in(2:end) - yp(1:end-1);


Theta    =  atan2(deltay,deltax);
Xhat1    =  deltax.*cos(Theta) + deltay.*sin(Theta);
Xhat1_in =  deltax_across1.*cos(Theta) + deltay_across1.*sin(Theta);
Yhat1_in = -deltax_across1.*sin(Theta) + deltay_across1.*cos(Theta);
Xhat0_in =  deltax_across.* cos(Theta) + deltay_across.* sin(Theta);
Yhat0_in = -deltax_across.* sin(Theta) + deltay_across.* cos(Theta);
m1       =  Yhat0_in./Xhat0_in;
m2       =  Yhat1_in./(Xhat1_in - Xhat1);
m3       =  (Yhat1_in - Yhat0_in)./(Xhat1_in - Xhat0_in);
c3       =  Yhat0_in - m3.*Xhat0_in;


raw = zeros(npts,o);
for i = 1:npts
    
	% Defining window that we care about 
    % Don't want to do all these calculations on the whole image
	[~,jx]  = roundx(xp(i),x_im);
	[~,iy]  = roundx(yp(i),y_im);
	j1      = max(jx-(depthInsideImage+10),1);
	j2      = min(jx+(depthInsideImage+10),n);
	i1      = max(iy-(depthInsideImage+10),1);
	i2      = min(iy+(depthInsideImage+10),m);
	
	% Transforming coordinates on our little local rectangle - image
	theta  = Theta(i);
	Xhat   =  (X(i1:i2,j1:j2) - xp(i))*cos(theta) + (Y(i1:i2,j1:j2) - yp(i))*sin(theta);
	Yhat   = -(X(i1:i2,j1:j2) - xp(i))*sin(theta) + (Y(i1:i2,j1:j2) - yp(i))*cos(theta);
        
	u  = Yhat < m3(i)*Xhat+c3(i) & Yhat >= 0;
	v  = Xhat > (1/m1(i))*Yhat & Xhat < (1/m2(i))*Yhat + Xhat1(i);
	bw = u & v;	
    for j = 1:o
		I1       = I(i1:i2,j1:j2,j);
		I1       = I1(bw);
		raw(i,j) = mean(I1(:));
    end
end


%
% Smoothing. 
%
t = zeros(npts+1,o);
for j = 1:o
	I1      = [raw(1:end-1,j);raw(:,j);raw(2:end,j)];
	smth1   = smooth(I1,round(0.02*npts));            % PB - avging over 1% of points
	t(:,j)  = smth1(npts:2*npts);
end
raw = raw([1:end,1],:);
s   = linspace(-1,1,npts+1)';
t   = interp1(ss,t,s);
raw = interp1(ss,raw,s);
%plot(s,t)



% Other attempts at extracting domain
%{
%{
non-vectorized version

x = xp;
y = yp;
x_in = xp_in;
y_in = yp_in;

tic
% Defining image coordinates
x_im = (1:n);              y_im = (1:m)';
X    = repmat(x_im,m,1);   Y    = repmat(y_im,1,n);

%
% Building raw data
%
raw2 = zeros(ns,o);
for i = 1:ns
    
	% Defining window that we care about 
    % Don't want to do all these calculations on the whole image
	[jx,jx] = roundx(x(i),x_im);
	[iy,iy] = roundx(y(i),y_im);
	j1      = max(jx-(Yhatmax+10),1);
	j2      = min(jx+(Yhatmax+10),n);
	i1      = max(iy-(Yhatmax+10),1);
	i2      = min(iy+(Yhatmax+10),m);
	
	% Transforming coordinates on our little local rectangle - image
	theta  = atan2(y(i+1)-y(i),x(i+1)-x(i));
	Xhat   =  (X(i1:i2,j1:j2) - x(i))*cos(theta) + (Y(i1:i2,j1:j2) - y(i))*sin(theta);
	Yhat   = -(X(i1:i2,j1:j2) - x(i))*sin(theta) + (Y(i1:i2,j1:j2) - y(i))*cos(theta);
    
    % Transforming coordinates of our points!
	xhat1    =  (x(i+1) - x(i))*cos(theta) + (y(i+1) - y(i))*sin(theta); % xhat should be ~ds.
% 	yhat1    = -(x(i+1) - x(i))*sin(theta) + (y(i+1) - y(i))*cos(theta); % yhat should be zero.
    xhat1_in =  (x_in(i+1) - x(i))*cos(theta) + (y_in(i+1) - y(i))*sin(theta);
    xhat0_in =  (x_in(i)   - x(i))*cos(theta) + (y_in(i)   - y(i))*sin(theta);
    yhat1_in = -(x_in(i+1) - x(i))*sin(theta) + (y_in(i+1) - y(i))*cos(theta);
    yhat0_in = -(x_in(i)   - x(i))*sin(theta) + (y_in(i)   - y(i))*cos(theta);
    
    m1 = yhat0_in/xhat0_in;
    m2 = yhat1_in/(xhat1_in - xhat1);
    m3 = (yhat1_in - yhat0_in)/(xhat1_in - xhat0_in);
    c3 = yhat0_in - m3*xhat0_in;
    
	u  = Yhat < m3*Xhat+c3 & Yhat >= 0;
	v  = Xhat > (1/m1)*Yhat & Xhat < (1/m2)*Yhat + xhat1;
	bw = u & v;	
    for j = 1:o
		I1       = I(i1:i2,j1:j2,j);
		I1       = I1(bw);
		raw2(i,j) = mean(I1(:));
    end
end
toc
%}





% Another solution
%{
tic
xmid = [xp,xp_in];
ymid = [yp,yp_in];
D       = 55;
xp2    = xp - xmid + D;
yp2    = yp - ymid + D;
xp2_in = xp_in - xmid + D;
yp2_in = yp_in - ymid + D;
xmid = round(mean(xmid,2));
ymid = round(mean(ymid,2));

%D    = round(10*mean(abs(diff(xp))));


% check solution
%{
i = 1;
xall = xmid(i)-D:xmid(i)+D;
yall = ymid(i)-D:ymid(i)+D;
leftx  = repmat(xall(1),111,1);
rightx = repmat(xall(end),111,1);
upy = repmat(yall(1),111,1);
downy = repmat(yall(end),111,1);
figure
imshow(I,[])
hold on
plot(xp(i),yp(i),'*')
plot(xp_in(i),yp_in(i),'*')
plot(xmid(i),ymid(i),'*')
plot(leftx,yall,'.')
plot(rightx,yall,'.')
plot(xall,upy,'.')
plot(xall,downy,'.')
%}


raw = zeros(ns,o);
for i = 1:ns
    
	% Extracting window that we care about 
    I1   = I(ymid(i)-D:ymid(i)+D, xmid(i)-D:xmid(i)+D,:);
    mask = roipoly(I1,[xp2(i),xp2(i+1),xp2_in(i),xp2_in(i+1)],[yp2(i),yp2(i+1),yp2_in(i),yp2_in(i+1)]);
    I2   = I1.*uint16(mask);
    
    % check solutions
    %{
    figure
    imshow(I2(:,:,1),[])
    hold on
    plot(xp2(i),yp2(i),'*')
    %plot(xp2(i+1),yp2(i+1),'*')
    plot(xp2_in(i),yp2_in(i),'*')
    plot(xmid(i),ymid(i),'*')
    %plot(xp2_in(i+1),yp2_in(i+1),'*')
    %}
    
    for j = 1:o
        I3 = I2(:,:,j);
		raw(i,j) = mean(I3(:));
    end 
end
toc
%}


m_across = (yp_in - yp)./(xp_in - xp);
m_in = (yp_in(1:end) - [yp_in(end-1); yp_in(1:end-1)])./(xp_in(1:end) - [xp_in(end-1); xp_in(1:end-1)]);
m_out = (yp(1:end) - [yp(end-1); yp(1:end-1)])./(xp(1:end) - [xp(end-1); xp(1:end-1)]);

c_across = yp    - m_across.* xp;
c_in     = yp_in - m_in.* xp_in; 
c_out    = yp    - m_out.* xp; 

Xhat1 =  [xp(2:end);xp(1)] - xp(1:end) + [yp(2:end);yp(1)] - yp(1:end);
Yhat1 = -[xp(2:end);xp(1)] - xp(1:end) + [yp(2:end);yp(1)] - yp(1:end);
Xhat0 =  xp(1:end) - [xp(end-1); xp(1:end-1)] +  yp(1:end) - [yp(end-1); yp(1:end-1)];
Yhat0 = -xp(1:end) - [xp(end-1); xp(1:end-1)] +  yp(1:end) - [yp(end-1); yp(1:end-1)];

Xhat1_in =  [xp_in(2:end);xp_in(1)] - xp_in(1:end) +  [yp_in(2:end);yp_in(1)] - yp_in(1:end);
Yhat1_in = -[xp_in(2:end);xp_in(1)] - xp_in(1:end) + [yp_in(2:end);yp_in(1)] - yp_in(1:end);
Xhat0_in =  xp_in(1:end) - [xp_in(end-1); xp_in(1:end-1)] + yp_in(1:end) - [yp_in(end-1); yp_in(1:end-1)];
Yhat0_in = -xp_in(1:end) - [xp_in(end-1); xp_in(1:end-1)] + yp_in(1:end) - [yp_in(end-1); yp_in(1:end-1)];



by = 3;
iodd  = 1:by:ns;
ieven = 2:by:ns;

mask1 = roipoly(I(:,:,1),[xp(iodd),xp_in(iodd)],[yp(iodd),yp_in(iodd)]);

for i=1:ns
    [x2,y2,BW,xi2,yi2] = roipoly(I(:,:,1),[xp(i),xp(i+1),xp_in(i),xp_in(i+1)],[yp(i),yp(i+1),yp_in(i),yp_in(i+1)]);
end




tic
raw = zeros(ns,o);
for i=1:ns
    mask = roipoly(I(:,:,1),[xp(i),xp(i+1),xp_in(i),xp_in(i+1)],[yp(i),yp(i+1),yp_in(i),yp_in(i+1)]); 
    mask = I.*uint16(mask);
    for j=1:o
       raw(i,j) = mean(mask(:));
    end
end
toc
%}










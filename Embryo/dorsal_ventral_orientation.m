function theta = dorsal_ventral_orientation(Xp,Yp,Zp,xc,yc,len,phi)
% Function to find orientation
% This function finds the z-orientation of the embryo using
% the coordinates of (xp,yp) for different z-stacks

load('orientation.mat')
close all

% set controls
%
From        = 3;
To          = 9;
numpoints   = 200;
len         = length(From:To);

% get info
%
nSlices = size(Xp,2);
npts    = size(Xp,1);


% Here we find a center in XY plane that would act as a center for all
% zslices.
count = 1;
phi = 0;
xc = 0;
yc = 0;
l = inf;
h = 0;
for i =From:To
    xp = Xp(:,i); yp = Yp(:,i); 

    % Fit to an ellipse to rotate the image later
    [l0,h,xc,yc,phi] = ellipse_fit(xp(1:end-1),yp(1:end-1));
    if phi > pi/2
        phi = phi - pi;
    end
    phi = phi + phi;
    xc = xc + xc;
    yc = yc + yc;
    if l0<l
        l = l0;
    end
    h = h + h;
end
xc = xc/len;
yc = yc/len;
phi = phi/len;
xgrid = linspace(-l,l,numpoints);



% Prep Data
count = 1;
for i =From:To
    xp = Xp(:,i); yp = Yp(:,i); 
 
    %check by plotting
    %{
    % foci of ellipse:
    xf1 = xc + (sqrt(((.5*l)^2)-((.5*h1)^2)))*cos(phi);
    yf1 = yc + (sqrt(((.5*l)^2)-((.5*h1)^2)))*sin(phi);
    xf2 = xc - (sqrt(((.5*l)^2)-((.5*h1)^2)))*cos(phi);
    yf2 = yc - (sqrt(((.5*l)^2)-((.5*h1)^2)))*sin(phi);
    figure
    plot(xp,yp,'*k',xc,yc,'og',xf1,yf1,'og',xf2,yf2,'og')
    hold on
    m = ((yf2 - yf1)/(xf2 - xf1));
    xspan = linspace(min(xp),max(xp));
    y  = m*(xspan - xf1) + yf1;
    plot(xspan, y);
    %}
    
    
    % rotate points to make major axis parallel to x axis:
     xp = xp - xc;  yp = yp - yc; 
    [xp, yp] = rotatePoints(pi-phi,xp, yp);
    
    
    % check by plotting
    %{
    % rotating foci
    xf1 = xf1-xc;     xf2 = xf2-xc;     yf1 = yf1-yc;   yf2 = yf2-yc;
    xc  = 0;          yc  = 0;
    [xf1, yf1] = rotatePoints(pi-phi,xf1, yf1);
    [xf2, yf2] = rotatePoints(pi-phi,xf2, yf2);
    figure
    plot(xp,yp,'*k',xc,yc,'og',xf1,yf1,'og',xf2,yf2,'og')
    axis equal
    hold on
    m = ((yf2 - yf1)/(xf2 - xf1));
    xplot = linspace(1.1*min(xp),1.1*max(xp));
    y  = m*(xplot - xf1) + yf1;
    plot(xplot, y);
    %}
    
    
    xp      = xp(1:end-1);
    yp      = yp(1:end-1);
    k1      = find(yp>0);
    k2      = find(yp<0); 
    ygrid1(:,count)  = interp1(xp(k1),yp(k1),xgrid)';
    ygrid2(:,count)  = interp1(xp(k2),yp(k2),xgrid)'; 
    
    count = count+1;
   
end



% Assign 3D coordinate points
zgrid  = (Zp(1,From:To) - Zp(1,end))';
% xgrid  = xgrid'*scalings(1);
% ygrid1 = ygrid1.*scalings(2);
% ygrid2 = ygrid2.*scalings(2);

Yplot       = [ygrid1; ygrid2];
[zlen,xlen] = size(Yplot);
Zplot       = repelem(zgrid,zlen)';
Xplot       = repmat([xgrid;xgrid],xlen,1);
Yplot       = Yplot(:);



% Get Orientation by first fitting the cross-sections to a circle then
% fitting their centers to a parabola!
Zc = zeros(numpoints,1);    Yc = Zc;    Xc = Zc;
for i=10:numpoints
    ycirc = [ygrid1(i,:),ygrid2(i,:)]';
    npts  = length(ycirc);
    xcirc = repmat(xgrid(i),1,npts)';
    zcirc = [zgrid; zgrid];
    
    [Yc(i),Zc(i),R(i)] = circfit(ycirc, zcirc);
    Xc(i) = xgrid(i);
    %h           = plotcircle(xc,yc,R);   
    
    % check by plotting
    %
    figure
    plot(ycirc,zcirc,'*b'); hold on;
    axis equal
    plot(Yc(i),Zc(i),'*r')
    plotcircle(Yc(i),Zc(i),R(i))
    %}
end


% Plotting circle centers to a parabola
i0 = ~isnan(Xc) & ~isnan(Yc) & ~isnan(Zc);
Xc = Xc(i0);    Yc = Yc(i0);    Zc = Zc(i0);




% Fit the centers to a parabola
f = fit([Xc,Yc],Zc,'poly11');


% Make final plots!
figure
plot3(Xplot,Yplot,Zplot,'.')    % Borders
axis equal; hold on
% plot3(Xc,Yc,Zc,'*')             % Centers of cross-sections
plot(f)                         % Best-fit plane
set(gca,'XTick',[], 'YTick', [],'ZTick',[])
axis 'off'
xlim([-110,110]);
ylim([-110,110]);
zlim([-5,50]);


% FInd angle between the fit-plane and the XY plane
vec1 = [-f.p10,-f.p01,1];
vec2 = [0,0,1];
theta = acosd(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));

end

function h = plotcircle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit);
    %hold off
end
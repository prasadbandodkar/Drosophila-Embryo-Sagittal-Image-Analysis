function [iPoles,Rad] = find_poles(xp, yp)
% This functions first fits the periphery to an ellipse. It rotates the 
% image and finds the max and min x-points of the poles. Next the function
% determines the cruvature at these two extremes to determine the anterior
% pole and the posteior pole.
% Input - 
%       I - The image I; used only if yesplot is true
%       xp,yp - x and y coordinates of the periphery of the embryo
%       i - Slice number
%       varargin - yesplot, filename, outfolder - all of these are used
%       only when yesplot is enabled. 
 

% Fit to an ellipse and rotate to make x axis approximately parallel to 
% the AP axis
[~, ~, ~, ~, phi] = ellipse_fit(xp, yp);
if phi > pi/2
    phi = phi - pi;
end
[xp2 ,yp2] = rotatePoints(pi-phi,xp, yp); 



% pad points
% x    = [flipud(xp2); xp2; xp2];
% y    = [flipud(yp2); yp2; yp2];
x = [xp2; xp2(2:end-1); xp2];
y = [yp2; yp2(2:end-1); yp2];
lenx = length(xp2);

% find the two end points and extract points for curvature calculations
nCurve  = round(length(xp2)/25);
[~,i1]  = min(xp2);      i1 = i1+lenx;
[~,i2]  = max(xp2);      i2 = i2+lenx;
xPole1  = x(i1-nCurve:i1+nCurve);        xPole2 = x(i2-nCurve:i2+nCurve);
yPole1  = y(i1-nCurve:i1+nCurve);        yPole2 = y(i2-nCurve:i2+nCurve);


[xP1,rad1] = findCurvature(xPole1,yPole1);
[xP2,rad2] = findCurvature(xPole2,yPole2); 
[~,i1]     = min(abs(xP1 - xp2));   
[~,i2]     = min(abs(xP2 - xp2));
if rad1 < rad2
    iPoles = [i1,i2]; 
    Rad    = [rad1,rad2];
else
    iPoles = [i2,i1];
    Rad    = [rad2,rad1];
end   

% check by plotting
%{
figure
plot(xp2,yp2)
axis equal; hold on
plot(xPole1,yPole1,'LineWidth',5,'Color','r')
plot(xPole2,yPole2,'LineWidth',5,'Color','b')
set(gca,'XTick',[], 'YTick', [])
plot(xp2(iPoles(1)),yp2(iPoles(1)),'*r','MarkerSize',28)
plot(xp2(iPoles(2)),yp2(iPoles(2)),'*b','MarkerSize',28)
axis 'off'
%}

end



%__________________sub-function to find curvature__________________________
function [xend,Rmin] = findCurvature(xpend,ypend)
    y       = linspace(min(ypend),max(ypend),100);
    p       = polyfit(ypend,xpend,2);
    x       = polyval(p,y);
    %curvter = abs(2*p(1)./(1+(2*p(1)*y + p(2)).^1.5));
    %radCurv = 1./curvter;
    radCurv = abs((1+(2*p(1)*y + p(2)).^2).^1.5./(2*p(1)));
      
    
    [Rmin,iend]   = min(radCurv);
    xend = x(iend);
end


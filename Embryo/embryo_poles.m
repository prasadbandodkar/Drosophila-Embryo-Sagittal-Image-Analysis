function [Poles,gof] = embryo_poles(Xp,Yp,Zp,phi)

Xp = Xp(:);
Yp = Yp(:);
Zp = Zp(:);

% rotate points so that the x axis is approximately parallel to the AP axis
%
phi     = deg2rad(phi);
[Xp,Yp] = rotatePoints(-phi, Xp, Yp);     % rotate anti-clockwise


% pick coordinates to fit surface 
%
percentpts  = 0.15;
xmin        = min(Xp);     
xmax        = max(Xp);
x15         = round(percentpts*(xmax-xmin));         % 15% points from the pole
v1          = Xp < xmin + x15;        
v2          = Xp > xmax - x15;
yend1       = Yp(v1);    xend1 = Xp(v1);    zend1 = Zp(v1);
yend2       = Yp(v2);    xend2 = Xp(v2);    zend2 = Zp(v2);
ft          = fittype( @(a,b,c,d,e,x,y) a*x.^2 + b*x + c*y.^2 + d*y + e, ...
                'independent', {'x', 'y'},'dependent', 'z');
    

% side: 1
%
[sf,gof1] = fit([zend1, yend1],xend1,ft);
zp1       = -sf.b/(2*sf.a);
yp1       = -sf.d/(2*sf.c);
xp1       = sf(zp1,yp1);
K1        = curvature3D(sf,[zp1,yp1,xp1]);

% check by plotting
%{
figure
plot(sf,[zend1, yend1],xend1)
hold on
plot3(zp1,yp1,xp1,'*r')
axis equal
%}


% side: 2
%
[sf,gof2] = fit([zend2, yend2],xend2,ft);
zp2      = -sf.b/(2*sf.a);
yp2      = -sf.d/(2*sf.c);
xp2      = sf(zp2,yp2);
K2       = curvature3D(sf,[zp2,yp2,xp2]);

% check by plotting
%{
figure
plot(sf,[zend2, yend2],xend2)
hold on
plot3(zp2,yp2,xp2,'*r')
axis equal
set(gca,'XTick',[], 'YTick', [],'ZTick',[])
axis 'off'
%}

% check by plotting: Full embryo
%{
figure
plot3(Zp,Yp,Xp,'*')
axis equal
hold on
plot3(zp2,yp2,xp2,'*r')
plot3(zp1,yp1,xp1,'*r')
set(gca,'XTick',[], 'YTick', [])
axis 'off'
%}



% check for anterior pole using curvature
%
if K1 > K2
   Poles = [[xp1,yp1,zp1];[xp2,yp2,zp2]];   % side 1 is anterior pole
   gof   = [gof1;gof2]; 
else
   Poles = [[xp2,yp2,zp2];[xp1,yp1,zp1]];   % side 2 is anterior pole
   gof   = [gof1;gof2];
end


% rotate poles back by phi
%
[Poles(:,1),Poles(:,2)] = rotatePoints(phi, Poles(:,1), Poles(:,2));



%________________subfuntion to calculate curvature_________________________
function K = curvature3D(sf,pt)
    c = coeffvalues(sf);
    xval = pt(1);   yval = pt(2);
    P = 2*c(1)*xval + c(2);
    Q = 2*c(3)*yval + c(4);
    R = 2*c(1);
    S = 0;
    T = 2*c(3);
    
    K = abs((R*T - S^2)/(1 + P^2 + Q^2));        % gaussian curvature 
end



%ft = fittype( @(a,b,x0,y0,z0,x,y) a*(x-x0).^2+ b*(y-y0).^2 + z0, ...
% 'independent', {'x', 'y'},'dependent', 'z');

% check for anterior pole using estimates by findPoles function
%{
AntPoleEst = [Xp{6}(IPoles(6,1))*scalings(1), Yp{6}(IPoles(6,1))*scalings(2)];
PosPoleEst = [Xp{6}(IPoles(6,2))*scalings(1), Yp{6}(IPoles(6,2))*scalings(2)];
[AntPoleEst(1) ,AntPoleEst(2)] = rotatePoints(pi-phi,AntPoleEst(1), AntPoleEst(2)); 
[PosPoleEst(1) ,PosPoleEst(2)] = rotatePoints(pi-phi,PosPoleEst(1), PosPoleEst(2)); 

dist1 = norm([xp1,yp1] - AntPoleEst);
dist2 = norm([xp1,yp1] - PosPoleEst);

if dist1 < dist2
   Poles = [[xp1,yp1,zp1];[xp2,yp2,zp2]];   % side 1 is anterior pole
   side1 = true;
else
   Poles = [[xp2,yp2,zp2];[xp1,yp1,zp1]];   % side 2 is anterior pole 
   side1 = false;
end
%}
















end
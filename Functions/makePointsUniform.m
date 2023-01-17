function [x,y,ss] = makePointsUniform(xp,yp,ns)

%
% Perimeter of the embryo
%
dxp = diff(xp);   dyp = diff(yp);     ds = sqrt(dxp.^2 + dyp.^2);
perim = sum(ds);
sp = [0;cumsum(ds)]; % cumulative pseudo arclength for original perim



%
% Changing the density of points around the embryo
%
if isodd(ns)
	ns = ns + 1;
end
s = perim/ns*ones(ns,1);
s = [0;cumsum(s)]; % cumulative pseudo arclength for new perim

npts = length(xp);
if length(unique(xp)) < npts
   i0 = 1e-3 + rand(npts,1)*(1e-2 - 1e-3);      % add a small amount of noise to ensure unique points
   xp = xp + i0;
end
if length(unique(yp)) < npts
   i0 = 1e-3 + rand(npts,1)*(1e-2 - 1e-3);      % add a small amount of noise to ensure unique points
   yp = yp + i0;
end
if length(unique(sp)) < npts
   i0 = 1e-3 + rand(npts,1)*(1e-2 - 1e-3);      % add a small amount of noise to ensure unique points
   sp = sp + i0;
   sp(1) = 0;
end


x = interp1(sp,xp,s);
y = interp1(sp,yp,s);
x(end) = x(1);
y(end) = y(1);



dx = diff(x);       dy = diff(y);       ds = sqrt(dx.^2 + dy.^2);
perim2 = sum(ds);   ss = 2*[0;cumsum(ds)]/perim2 - 1;

% older method
%{
x = zeros(ns+1,1);
y = zeros(ns+1,1);
x(1) = xp(1); 
y(1) = yp(1);
Sr = 0; % running total of pseudoarclength.
i = 1;
for j  = 1:ns-1
	D  = sp(j+1) - Sr;
	D2 = D + sqrt((x(j)-xp(i))^2 + (y(j)-yp(i))^2);
	
	if D2 >= ds(i) && i < length(ds)
		D2 = D - sqrt((xp(i+1)-x(j))^2 + (yp(i+1)-y(j))^2);
		i  = i + 1;
	end
	d      = D2/ds(i);
	x(j+1) = d*dxp(i) + xp(i);
	y(j+1) = d*dyp(i) + yp(i);
	Sr     = Sr + sqrt((x(j+1) - x(j))^2 + (y(j+1) - y(j))^2);
end
x(end) = x(1); y(end) = y(1);

dx = diff(x);       dy = diff(y);       ds = sqrt(dx.^2 + dy.^2);
perim2 = sum(ds);   ss = 2*[0;cumsum(ds)]/perim2 - 1;
%}











function [data,IM] = analyze_embryo(IM,data)

% load('embryo.mat')
% close all

% extract embryo level data
%
Xp       = data.border.Xp;
Yp       = data.border.Yp;
scalings = data.metadata.scalings;
nSlices  = data.nSlices;
npts     = data.npts;


i0 = 1:nSlices;

% move from image coordinates to real coordinates
%
Xp = Xp.*scalings(1);
Yp = Yp.*scalings(2);
Zp = repmat(fliplr(1:nSlices).*scalings(3),npts,1);
Xp = Xp(1:npts,i0);
Yp = Yp(1:npts,i0);
Zp = Zp(1:npts,i0);


% check by plotting
%{
figure
plot3(Xp2,Yp2,Zp,'*')
axis equal
%}


% Task 1: find positions of the anterior and posterior poles in real coordinates
%
APphi       = data.metadata.APphi;
[Poles,gof] = embryo_poles(Xp,Yp,Zp,APphi);


% Task 2: find dorsal-ventral orientation of the embryo
%
theta = [];
APxc  = data.metadata.APxcenter*scalings(1);
APyc  = data.metadata.APycenter*scalings(2);
APlen = data.metadata.APlen*scalings(1);        % assuming scaling in x and y direction is equal
% theta = dorsal_ventral_orientation(Xp,Yp,Zp,APxc,APyc,APlen,APphi);


% store in data
%
data.Poles      = Poles;
data.DVAngle    = theta;
data.metadata.PolesGof = gof;


end
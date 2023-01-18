function [data,IM] = analyze_images(IM,data)
% First we find the periphery of the embryo and detect the intensity of
% each channel as you go around the periphery. We go through a loop that 
% processed each slice individually.
%


nSlices      = data.nSlices;
npts         = data.npts;
depthInImage = data.metadata.depthInImage;
channelnames = data.channelnames;
nChannels    = length(channelnames);
scalings     = data.metadata.scalings;



% run analysis over all images in the z-stack
%
IPoles            = zeros(nSlices,2);      
RadiusOfCurvature = IPoles;
Xp                = zeros(npts+1,nSlices);
Yp                = Xp;
for i=1:nSlices
    I  = squeeze(IM(:,:,:,i)); 
        
    % Find borders
    [xp,yp]     = border_finderI(I,npts);
    Xp(:,i)     = xp;
    Yp(:,i)     = yp;

    % Find poles and gene expression
    [iPoles,rad]  = find_poles(xp,yp); 
    [t,raw,s]     = domain_measure(I,xp,yp,depthInImage,npts);           % Find values of intensities for all channels
    
    % Store Values
    IPoles(i,:)             = iPoles;
    RadiusOfCurvature(i,:)  = rad;

    for j=1:nChannels
        name = channelnames(j);
        data.(name).T(:,i) = t(:,j);
        data.(name).meta.Raw(:,i) = raw(:,j);
    end
end


% find angle to rotate images so that ap axis is approximately parallel
% to the x axis
%
I                = max(IM,[],3);
I                = squeeze(max(I,[],4));
[xp,yp]          = border_finderI(I,npts);
[len,~,xc,yc,phi] = ellipse_fit(xp, yp);
if phi > pi/2
    phi = phi - pi;
end
phi = rad2deg(phi);


% store data
% 
ArcImage                        = sum(sqrt(diff(Xp).^2 + diff(Yp).^2));
ArcEmbryo                       = ArcImage/scalings(1);
data.border.Xp                  = Xp;
data.border.Yp                  = Yp;
data.border.IPoles              = IPoles;
data.border.ArcImage            = ArcImage;
data.border.ArcEmrbyo           = ArcEmbryo;
data.metadata.RadiusOfCurvature = RadiusOfCurvature;
data.metadata.APphi             = phi;
data.metadata.APlen             = len;
data.metadata.APxcenter         = xc;
data.metadata.APycenter         = yc;

























end
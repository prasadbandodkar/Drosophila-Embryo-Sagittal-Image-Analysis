function [IM,data,metadata] = image_pre_processing(IM,metadata,options)
% This function first extracts the region of interest as specified in the
% metadata and then performs background subtraction


% read user-provided options
%
depthInEmbryo = options.depthInEmbryo;
path_data     = options.path_data;
npts          = options.npts;

% check channel information and remove images that will not be analyzed
%
fileloc = metadata.fileloc;
[filename,genotype,channelinf] = parse_filename(fileloc);
nChannels0  = metadata.lsminf1.DimensionChannels;
nChannels   = length(channelinf.channelnames);
if nChannels0 ~= nChannels
    warning(['The number of channels parsed from the filename did not' ...
    ' match that found in the metadata. Continuting analysis with' ...
    ' user provided information'])
end
i0 = channelinf.channelIDs ~= 5;
nChannels = sum(i0);
IM(:,:,~i0,:) = [];             % Remove channel 5 - DIC image
[H,W,nChannelsI,nSlices] = size(IM);  
if nChannelsI ~= nChannels
    warning(['The number of channels found in the image data did not' ...
    ' match that found in the metadata. Continuting analysis with' ...
    ' image data'])
end
nChannels = nChannelsI;


%
% extracting our region of interest:
% There will be a rectangle that contains the embryo, outside of which all 
% pixels will be zero.
%
lsminf2 = metadata.lsminf2;
if lsminf2.ScanInfo.USE_ROIS && ~lsminf2.ScanInfo.USE_REDUCED_MEMORY_ROIS
	I     = sum(sum(IM,3),4);
	maxbw = logical(max(I));    j = find(maxbw);  j1 = j(1); j2 = j(end);
	bwj   = logical(I(:,j1));   i = find(bwj);    i1 = i(1); i2 = i(end);
	IM    = IM(i1:i2,j1:j2,:,:);
end


% background subtraction on the images
%
[IM,bg,bgsig]   = subtract_background_image(IM,nChannels); 


% get scaling information
%
scalings     = 1e6*[metadata.lsminf1.VoxelSizeX,metadata.lsminf1.VoxelSizeY,metadata.lsminf1.VoxelSizeZ];   % microns/pixel
depthInImage = round(depthInEmbryo/scalings(1));         % about 30 pixels
rho          = round(scalings(3)/scalings(1));           % hope we don't have to round much.



% store metadata
%
metadata.bg            = bg;
metadata.std_bg        = bgsig;
metadata.scalings      = scalings;
metadata.rho           = rho;
metadata.depthinEmbryo = depthInEmbryo;
metadata.depthInImage  = depthInImage;


% store data
%
data.filename     = filename;
data.genotype     = genotype;
data.fileloc      = fileloc;
data.channelnames = channelinf.channelnames(i0);
data.channeltypes = channelinf.channeltypes(i0);
data.channelIDs   = channelinf.channelIDs(i0);
data.H            = H;
data.W            = W;
data.nSlices      = nSlices;
data.npts         = npts;
data.metadata     = metadata;

end
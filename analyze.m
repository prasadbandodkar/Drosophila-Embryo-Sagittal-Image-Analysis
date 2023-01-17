function data = analyze(filelocation,path_data,config,varargin)


% unpack varargin
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	channelID = varargin{iArg};
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	channelnames = varargin{iArg};
end,iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	genotype = varargin{iArg};
end
if nArg == 0
    c         = strsplit(filelocation,{filesep,'.'});
    filename  = char(c(end-1));
    genotype  = char(c(end-2));
    c         = strsplit(filename,' ');  
    [channelID, channelnames]    = matchChannels(c(end)); 
else
    c         = strsplit(filelocation,{filesep,'.'});
    filename  = char(c(end-1));
end


% extract configuration
%
if exist('configuration','var')
    depthInEmbryo   = config.depthInsideEmbryo;
    yesplot         = config.yesplot;
    npts            = config.npts;
else
    depthInEmbryo   = 18.36;      % in microns
    yesplot         = false;
    npts            = 2000;
end


% specify channel types
%
i0               = channelID == 5;
channelnames(i0) = [];
channelID(i0)    = [];
nChannels        = length(channelID);
channeltype      = strings(1,nChannels);
sides            = struct('sPeaks',[],'AntBorder',[],'PosBorder',[],'Width',[],'Gof',[], ...
                          'dsPeaks',[],'dAntBorder',[],'dPosBorder',[],'dWidth',[]);
meta             = struct('Raw',[],'side1',sides,'side2',sides); 
mRNA             = struct('channeltype', [],'channelID',[],'nPeaks',[],'sPeaks',[], ...
                          'AntBorder',[],'PosBorder',[],'Width',[],'T',[],'s',[],'meta',meta);

% Note: For other channels, the structure can change accordingly. The
% fitting step currently only works for the mRNA channel. 
for i = 1:nChannels
	switch channelID(i)
		case 1
			channeltype(i) = "nuclei";
            name = channelnames(i);
            chdata.(name).channeltype = channeltype(i);
            chdata.(name).channelID   = 1;
            chdata.(name) = mRNA;
		case 2
			channeltype(i) = "non-nuclear protein or mRNA";
            name = channelnames(i);
            chdata.(name).channeltype = channeltype(i);
            chdata.(name).channelID   = 2;
            chdata.(name) = mRNA;
		case 3
			channeltype(i) = "nuclear protein";
            name = channelnames(i);
            chdata.(name).channeltype = channeltype(i);
            chdata.(name).channelID   = 3;
            chdata.(name) = mRNA;
		case 4
			channeltype(i) = "intronic probe";
            name = channelnames(i);
            chdata.(name).channeltype = channeltype(i);
            chdata.(name).channelID   = 4;
            chdata.(name) = mRNA;
		case 5
			channeltype(i) = "N/A";
		otherwise
			channelID(i) = 5;
			channeltype(i) = "N/A";
	end
end
nuc_ch = find(channelID == 1);
if length(nuc_ch) > 1
	error('You can specify only one nuclear channel.')
end


% early remarks
%
data.filename               = filename;
data.filelocation           = filelocation;
data.channelnames           = channelnames;
data.channeltypes           = channeltype;
data.channels               = channelID;
data.genotype               = genotype;



% reading raw data (currently only supports reading lsm files)
% 
lsm                  = lsmRead2(filelocation);
lsminf1              = lsm.lsm;
lsminf2              = lsminfo(filelocation);
scalings             = 1e6*[lsminf1.VoxelSizeX,lsminf1.VoxelSizeY,lsminf1.VoxelSizeZ];   % microns/pixel
depthInImage         = round(depthInEmbryo/scalings(1));         % about 30 pixels
rho                  = round(scalings(3)/scalings(1));           % hope we don't have to round much.
IM                   = lsm.data;   
IM(:,:,i0,:)         = [];                                       % Remove channel 5 - DIC image
[H,W,nChannelsI,nSlices] = size(IM);                             % height, width, channels, zstacks
clear lsm     
if nChannels~=nChannelsI
    fprintf('#channels specified: %d\n#channels found: %d\n',nChannels,numchI)
    error('The number of channels specified by user does not match that found in the file')
end


% subtract background
%
[IM,bg,bgsig] = subtract_background_image(IM,nChannels); 


% extracting our region of interest.  
%
% There will be a rectangle that contains the embryo, outside of which all 
% pixels will be zero.
if lsminf2.ScanInfo.USE_ROIS && ~lsminf2.ScanInfo.USE_REDUCED_MEMORY_ROIS
	I     = sum(sum(IM,3),4);
	maxbw = logical(max(I));    j = find(maxbw);  j1 = j(1); j2 = j(end);
	bwj   = logical(I(:,j1));   i = find(bwj);    i1 = i(1); i2 = i(end);
	IM    = IM(i1:i2,j1:j2,:,:);
end


% middle remarks:
%
data.H                = H; 
data.W                = W; 
data.nSlices          = nSlices;
data.config           = config;
metadata.lsminf1      = lsminf1;
metadata.lsminf2      = lsminf2;
metadata.scalings     = scalings;
metadata.rho          = rho;
metadata.bg           = bg;
metadata.std_bg       = bgsig;
metadata.depthInImage = depthInImage;






% Analysis
%
% First we find the periphery of the embryo and detect the intensity of
% each channel as you go around the periphery. We go through a loop that 
% processed each slice individually.
%
% initialize variables
%
IPoles            = zeros(nSlices,2);      
RadiusOfCurvature = IPoles;
Xp                = zeros(npts+1,nSlices);
Yp                = Xp;

% loop over zstack
%
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
        chdata.(name).T(:,i) = t(:,j);
        chdata.(name).meta.Raw(:,i) = raw(:,j);
    end
end


% find angle to rotate images so that ap axis is approximately parallel
% to the x axis
I                = max(IM,[],3);
I                = squeeze(max(I,[],4));
[xp,yp]          = border_finderI(I,npts);
[~,~,~,~,phirot] = ellipse_fit(xp, yp);
if phirot > pi/2
    phirot = phirot - pi;
end
phirot = rad2deg(phirot);
data.APangle = phirot;


if yesplot
    plot_raw_expression(IM,channelnames,path_data)
    plot_boundaries_poles(IM,Xp,Yp,phirot,IPoles,path_data);
end



% final remarks
%
ArcImage              = sum(sqrt(diff(Xp).^2 + diff(Yp).^2));
ArcEmbryo             = ArcImage/scalings(1);
data.border.Xp        = Xp;
data.border.Yp        = Yp;
data.border.IPoles    = IPoles;
data.border.ArcImage  = ArcImage;
data.border.ArcEmrbyo = ArcEmbryo;
metadata.Poles.RadiusOfCurvature   = RadiusOfCurvature; 
for i=1:nChannels
    name = channelnames(i);
    data.(name) = chdata.(name);
end
data.metadata = metadata;
































end
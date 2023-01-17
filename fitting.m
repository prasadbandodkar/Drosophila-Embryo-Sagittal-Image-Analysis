function data = fitting(data,varargin)
% 
% close all
% load('fitting.mat')

% unpack varargin
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	path_data = varargin{iArg}; else
	path_data = './';
end, iArg = iArg + 1; 
if nArg >= iArg && ~isempty(varargin{iArg}) 
	bkgrndthresh = varargin{iArg}; else
	bkgrndthresh = 0.15;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	peakthresh = varargin{iArg}; else
	peakthresh = 0.1;
end %, iArg = iArg + 1;


% extract info from data
%
filename     = data.filename; 
genotype     = data.genotype;
channelnames = data.channelnames;	
channelID    = data.channels;
nSlices      = data.nSlices;
IPoles       = data.border.IPoles;


% select data for fitting
%
vCh          = find(channelID == 2);      % "2" is the designation for mRNA
nChannels    = length(vCh);             % Total number of channels for fitting
channelnames = channelnames(vCh);


% Fitting
%
AllGeneInfo = load("AllGeneInfo.mat");
AllGeneInfo  = AllGeneInfo.AllGeneInfo;
for i=1:nChannels
    channelname = char(channelnames(i));
   
    % get smoothened data and background info
    %
    T      = data.(channelname).T;
    rbr    = max(AllGeneInfo.(channelname).rbr);
    nPeaks = AllGeneInfo.(channelname).nPeaks;
    s0     = AllGeneInfo.(channelname).xPeakLocation;

    [fitdata, meta1, meta2] = fit_channel(T,rbr,channelname,nPeaks,s0,IPoles);

    data.(channelname).sPeaks    = fitdata.sPeaks;
    data.(channelname).AntBorder = fitdata.AntBorder;
    data.(channelname).PosBorder = fitdata.PosBorder;
    data.(channelname).Width     = fitdata.Width;
    

    % add meta1 and meta2 to data
    %
    fn = fieldnames(meta1.side1);
    for j=1:length(fn)
        data.(channelname).meta.side1.(fn{j}) = meta1.side1.(fn{j});
        data.(channelname).meta.side2.(fn{j}) = meta1.side2.(fn{j});
    end
    fn = fieldnames(meta2.side1);
    for j=1:length(fn)
        data.metadata.side1.(fn{j}) = meta2.side1.(fn{j});
        data.metadata.side2.(fn{j}) = meta2.side2.(fn{j});
    end
    
end


data.yesfit = true;



















































end
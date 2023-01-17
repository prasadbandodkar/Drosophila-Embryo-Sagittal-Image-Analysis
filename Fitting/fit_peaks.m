function [data,meta] = fit_peaks(t,channelname,nPeaks,H,I,bcoor,s0)


if ~exist('H','var')
    [H,I,bcoor] = find_peak_locations(t,nPeaks,s0);
end


%
% Looping over all peaks
%
hh     = 0.5;
ns     = length(t);
s      = linspace(0,1,ns)';
tFit   = 0;
for i=1:nPeaks
    
    % extract peak data
    ti = zeros(ns,1); 
    ti(bcoor(i,1):bcoor(i,2)) = t(bcoor(i,1):bcoor(i,2));

    % fitting
    if nPeaks > 1
        genename = [channelname,'_',num2str(i)];
    else
        genename = channelname;
    end
    ii = num2str(i);
    params = ['A',ii,'a',',delt',ii,'a',',x0',ii,'a'];
    fittypestr      = ['genefit(''',genename,''',x,B,',params,')'];
    f               = fittype(fittypestr);
    [cfun,cvals,cint,cint68,rsquare] = fitelephant(ti,s,H(i),I(i),f);


    % Getting the borders of the canonical gene representation
    % Note -  "A": anterior, "P": posterior
    [xA,xP,x_offset] = canonicalgeneborders(genename,hh);
    if isnan(xA), xA = 0; end
    if isnan(xP), xP = 1; end

    % Unpacking cvals and cints (goes in metadata)
    A(i)       = cvals(1);                  
    B(i)       = cvals(2);           
    delt       = cvals(3);                       
    x0         = cvals(4);            
    dcint      = 0.5*diff(cint68);    
    ddelt      = dcint(3);          
    dx0        = dcint(4);

    % Calculating Anterior and Posterior boundaries: 
    % Calculating borders of this particlar gene expression peak
    data.sPeaks(i)      = x0;
    data.AntBorder(i)   = (xA-x_offset)*delt + x0; 
    data.PosBorder(i)   = (xP-x_offset)*delt + x0;
    data.Width(i)       = data.PosBorder(i) - data.AntBorder(i);
    data.dsPeaks(i)     = dx0;
    data.dAntBorder(i)  = sqrt((xA-x_offset)^2*ddelt.^2 + dx0.^2);
    data.dPosBorder(i)  = sqrt((xP-x_offset)^2*ddelt.^2 + dx0.^2);
    data.dWidth(i)      = data.dPosBorder(i) - data.dAntBorder(i); 


    % Storing important fit-related variables      
    meta.Delt(i)        = delt;        
    meta.Rsquare(i)     = rsquare; 
    meta.Cint{i}        = cint;        
    meta.Cint68{i}      = cint68;      

    % Constructing fit profile
    tFit = tFit + cfun(s);     
end


% store the final tFit into the metadata.
%
meta.tFit = tFit;


% Scoring
%
% score is calculated using data from the anterior border of the first peak
% to the posterior border of the last peak
%
[~, iLeft]  = min(abs(s - data.AntBorder(1)));
[~, iRight] = min(abs(s - data.PosBorder(nPeaks)));
i0    = iLeft:iRight;
score = sqrt(sum((tFit(i0) - t(i0)).^2))/length(t(i0));




















































end
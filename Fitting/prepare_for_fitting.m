function T = prepare_for_fitting(T,channelname)


channelnames = data.channelnames;	
channelID    = data.channelID;

% constructing arrays for selection
%
k_mrna       = find(channelID == 2);      % "2" is the designation for mRNA
kCh          = sort(k_mrna);
nChannels    = length(kCh);             % Total number of channels for fitting


genenames = channelnames(k_mrna);
for j=1:nChannels
    genename = char(genenames(j));
    v        = contains(GeneInfoMat,genename);
    nPeaks   = AllGeneInfo(v).nPeaks;
    s0              = zeros(n_peaks,1);
    sAnt            = s0;
    sPos            = s0;
    use_x0          = true(n_peaks,1);
    rbr0            = s0;
    

end




end
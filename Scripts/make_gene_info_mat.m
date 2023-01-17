% This script is used to generate AllGeneInfo.mat. 
% This mat file is then used in fit_peaks!

clear
close all
clc

% Add gene names and n_peaks for each gene here
GeneNames = {'Kr';'Eve'};
nPeaks    = [1,7];

for i=1:length(GeneNames)
    genename = char(GeneNames(i));
    AllGeneInfo.(genename).nPeaks    = nPeaks(i);
    for j=1:nPeaks(i)
        if nPeaks(i) > 1
           loadname = [genename,'_',num2str(j),'avg.mat'];
        else
            loadname = [genename,'avg.mat'];
        end
        load(loadname,'s','t','rbr','s_peak','sD','sV')
        AllGeneInfo.(genename).(['s',num2str(j)]) = s;
        AllGeneInfo.(genename).(['t',num2str(j)]) = t;
        AllGeneInfo.(genename).rbr(j)    = rbr;
        AllGeneInfo.(genename).xPeakLocation(j) = s_peak;
        AllGeneInfo.(genename).xAnteriorBorder(j) = sV;
        AllGeneInfo.(genename).xPosteriorBorder(j) = sD;      
    end 
    AllGeneInfo.(genename) = orderfields(AllGeneInfo.(genename));
end

save("Mats/AllGeneInfo.mat","AllGeneInfo")
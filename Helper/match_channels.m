function [channel, channelnames] = match_channels(name)
% This function finds the names of genes and macthes them to the
% corresponding channel

    load("Mats/MatchChannels.mat","geneinfo","filenames")
    c               = strcmp(name,filenames);
    if all(c==false,'all')
        channel = [];
        channelnames = [];
        return
    end
    channelnames    = strsplit(filenames(c,2)," ");
    channelnames    = [channelnames, "DAPI", "DIC"];
    channel         = zeros(1,length(channelnames));
    pos             = channel;
    
    count = 1;
    for i=1:length(channelnames)
       c = strcmp(channelnames(i),geneinfo.genenames);
       if any(c)
           channel(count) = geneinfo.channels(c);
           pos(count)     = geneinfo.positions(c);
           count          = count + 1;
       end
    end
    
    [~,pos]      = sort(pos);
    channelnames = channelnames(pos);
    channel      = channel(pos);
    
end
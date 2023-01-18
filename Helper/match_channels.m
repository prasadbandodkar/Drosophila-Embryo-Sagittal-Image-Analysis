function [channelIDs, channelnames, channeltypes] = match_channels(name)
% This function finds the names of genes and macthes them to the
% corresponding channel. Note that the assignment here is made solely based
% on parsing the filename. The true number of channels is obtained from the
% metadata of the raw image file.

    load("Mats/MatchChannels.mat","geneinfo","filenames")

    % get channelnames from filename
    %
    c  = strcmp(name,filenames);
    if all(c==false,'all')
        channelIDs = [];
        channelnames = [];
        return
    end
    channelnames    = strsplit(filenames(c,2)," ");
    channelnames    = [channelnames, "DAPI", "DIC"];        % Add DAPI and DIC
    channelIDs      = zeros(1,length(channelnames));
    pos             = channelIDs;
    

    % get channelIDs from the matfile
    %
    count = 1;
    for i=1:length(channelnames)
       c = strcmp(channelnames(i),geneinfo.genenames);
       if any(c)
           channelIDs(count) = geneinfo.channels(c);
           pos(count)     = geneinfo.positions(c);
           count          = count + 1;
       end
    end
    [~,pos]      = sort(pos);
    channelnames = channelnames(pos);
    channelIDs   = channelIDs(pos);


    % get channeltypes from channelIDs
    %
    nChannels = length(channelIDs);
    for i = 1:nChannels
	    switch channelIDs(i)
		    case 1
			    channeltypes(i) = "nuclei";
		    case 2
			    channeltypes(i) = "non-nuclear protein or mRNA";
		    case 3
			    channeltypes(i) = "nuclear protein";
		    case 4
			    channeltypes(i) = "intronic probe";
		    case 5
			    channeltypes(i) = "N/A";
            otherwise
			    channeltypes(i) = "N/A";
	    end
    end
    
end
function [filename,genotype,channelinf] =  parse_filename(file)
% This function extracts filename and genotype by parsing the location of
% the raw image. It also creates directories to store the final output in. 

warning('off', 'MATLAB:MKDIR:DirectoryExists');


% get filename and genotype from filename
%
c        = strsplit(file,{filesep,'.'});
try
    filename = char(c(end-1));
    genotype = char(c(end-2));   
catch
    filename = file;
    genotype = '';
end


% get channel information 
%
c = strsplit(filename,' ');  
[channelIDs, channelnames, channeltypes] = match_channels(c(end));
channelinf.channelIDs   = channelIDs;
channelinf.channelnames = channelnames;
channelinf.channeltypes = channeltypes;

end
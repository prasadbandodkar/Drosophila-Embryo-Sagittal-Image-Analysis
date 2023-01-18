function [IM,metadata] = read_raw_image(file,metadata)

arguments
    file string
    metadata = []
end

% reading raw data (currently only supports reading lsm files)
% 
fmt = strsplit(file,'.');
fmt = fmt{end};

if strcmp(fmt,'lsm')
    lsm              = lsmRead2(file);
    IM               = lsm.data;      
    metadata.fileloc = file;
    metadata.lsminf1 = lsm.lsm;
    metadata.lsminf2 = lsminfo(file);
else
    error('File format is currently not supported')
end


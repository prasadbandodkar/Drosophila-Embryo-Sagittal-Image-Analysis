function Files = get_files(path_raw,path_data,mode,filetype)
    % If the mode is 'test', the function returns 50 randomly selected
    % filenames from path_raw. If the mode is 'production', the functions
    % returns the filenames that haven't been processed yet in path_data

    arguments
        path_raw
        path_data
        mode     = 'production'
        filetype = 'lsm'
    end

    
    % Get lsm file locations
    Files   = extractFileLocations(path_raw,filetype);
    
    
    % check mode
    if ~isempty(Files)
        if strcmp(mode,'test')
            rng(0)
            N       = 50;     % if in test mode, code runs N random lsm files
            nFiles  = length(Files);   
            i0      = 1 + round(rand(N,1).*(nFiles - 1));
            Files   = Files(i0);
        elseif strcmp(mode,'production')
            Files   = getUnprocessedFiles(Files,path_data);
        end
    end


end
function run_analysis(path_raw, path_data, config, opts)
% This functions runs analysis over all files in path_raw and stores data
% in path_data using the configuration stored in configuration. The mode
% can be either 'test' or 'production'. When it's 'test', it runs the code
% over 50 randomly selected files. When it's 'production' it runs over all
% files. Default mode is 'production'.

arguments 
    path_raw 
    path_data 
    config                    = []
    opts.mode                 = 'production'
    opts.overwriteDatabase    = false;
    opts.yesplot              = true;
end

fprintf('Analyzing files in %s mode\n\n',opts.mode)


warning('off', 'MATLAB:MKDIR:DirectoryExists');
if ~isfolder(path_data)
    mkdir(path_data)
end


% get files
%
Files = get_files(path_raw, mode=opts.mode);


% change this in the future. Runs only those files that contain Kr and Eve.
Files   = Files(contains(Files,'BcdKrEve'));



%
% Loop over all files
%
nFiles    = length(Files); 
for i = 2:nFiles
    
    tic    

    % get filename and genotype. Make folders, if necessary.
    %
    file     = Files{i};
    c        = strsplit(file,{filesep,'.'});
    filename = char(c(end-1));
    genotype = char(c(end-2));   
    path_gen = [path_data,genotype];
    if ~isfolder(path_gen)
        mkdir(path_gen)
    end
    path_datafile = [path_gen,'/',filename];
    mkdir(path_datafile)
    

    % check if the filename exists on FlySection database

    
    % get channels and channelnames from the filename.
    c = strsplit(filename,' ');  
    [channelID, channelnames] = match_channels(c(end));

    
    % run analysis
    fprintf('Running:%d/%d %s\n',i,nFiles, filename)
     try    
        data = analyze(file,path_datafile,config,channelID,channelnames,genotype);  
        data = fitting(data,path_datafile); 

        % upload data to database
        if opts.overwriteDatabase
            upload2Firebase(data);
        end
        
        % save file. This fn is helpful when running parallel for loop
        parsave([path_datafile,'/',filename],data)
       
    catch ME
        disp(['Error message: ',ME.message])
        fprintf('Error location : %s in line: %d \n', ME.stack(1).name,ME.stack(1).line)
    end
   
    toc 
    
    close all
    fprintf('\n')
end


end


























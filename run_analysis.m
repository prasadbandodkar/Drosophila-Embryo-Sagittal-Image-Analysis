function run_analysis(path_raw, options)
% This functions runs analysis over all files in path_raw and stores data
% in path_data using the configuration stored in configuration. The mode
% can be either 'test' or 'production'. When it's 'test', it runs the code
% over 50 randomly selected files. When it's 'production' it runs over all
% files. Default mode is 'production'.

arguments 
    path_raw 
    options.mode                 = 'production'
    options.overwriteDatabase    = false;
    options.yesplot              = true;
    options.depthInEmbryo        = 30;
    options.npts                 = 2000;
    options.path_data            = path_raw;
end

fprintf('Analyzing files in %s mode\n\n',options.mode)

path_data = options.path_data;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
if ~isfolder(path_data)
    mkdir(path_data)
end


% get files
%
Files = get_files(path_raw, mode=options.mode);


% change this in the future. Runs only those files that contain Kr and Eve.
Files   = Files(contains(Files,'BcdKrEve'));



%
% Loop over all files
%
nFiles    = length(Files); 
for i = 1:nFiles
    file     = Files{i};
    fprintf('Running:%d/%d %s\n',i,nFiles, file)

    tic    
    try 
        % processing image
        %
        [IM,metadata] = read_raw_image(file);
        [IM,data]     = image_pre_processing(IM,metadata,options);
        [data,IM]     = analyze_images(IM,data);

        % processing embryo
        %
        [data,IM] = analyze_embryo(IM,data); 

        % fitting
        %
        data = fitting(data,path_datafile); 

        % upload data to database
        if options.overwriteDatabase
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


























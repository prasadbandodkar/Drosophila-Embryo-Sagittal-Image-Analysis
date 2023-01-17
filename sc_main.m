% script_main
% This is the main script to analyze the DGRP image dataset. 
clear
clc
close all

% addpath(genpath('./'))
 

% path configuration
path_raw          = '/Volumes/Extreme/Projects/DGRP_Project/Data_lsm/';
path_data         = '../Data/10/';


% options
yesplot = true;
mode    = 'test';


% Analysis configuration
config.depthInEmbryo = 30;                 % in microns. This was originally 18.36
config.npts          = 4000;


run_analysis(path_raw,path_data,config,mode=mode,yesplot=yesplot);






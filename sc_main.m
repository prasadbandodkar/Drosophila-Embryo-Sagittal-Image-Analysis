% script_main
% This is the main script to analyze the DGRP image dataset. 
clear
clc
close all

% addpath(genpath('./'))
 

% path configuration
path_raw          = '/Volumes/Extreme/Projects/DGRP_Project/Data_lsm/';
path_data         = '../Data/10/';



run_analysis(path_raw,path_data=path_data,yesplot=true,depthInEmbryo=30,npts=2000,mode='test');






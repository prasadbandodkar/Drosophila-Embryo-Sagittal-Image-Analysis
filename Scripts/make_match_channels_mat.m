close all
clear
clc

genenames   = ["Bcd","Eve","Gt","Hb","Kr","Kni","DAPI","DIC"];
channels    = [3,2,2,2,2,2,1,5];
positions   = [5,1,2,0,2,0,4,3];        % Only works for files that have Bcd, Kr, Eve



geneinfo = table(genenames',channels',positions','VariableNames',{'genenames','channels','positions'});
filenames = ["BcdGtEve","Bcd Gt Eve";
             "BcdKrEve", "Bcd Kr Eve"; 
             "HbGtKni", "Hb Gt Kni"];

save("Mats/MatchChannels.mat","geneinfo","filenames")
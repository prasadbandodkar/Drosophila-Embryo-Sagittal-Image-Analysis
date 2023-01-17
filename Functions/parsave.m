function parsave(savename,data)
% save files when running parallely
    save([savename,'.mat'],'data')
end
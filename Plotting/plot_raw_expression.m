function plot_raw_expression(IM,channelnames,path_data)

figure
nChannels = length(channelnames);
for i=1:nChannels
    I  = squeeze(IM(:,:,i,:));
%     I  = imrotate(I,phirot);
    montage(I,'DisplayRange',[])

    title(channelnames(i))
    saveas(gcf,[path_data,'/',char(channelnames(i)),'.tiff'])
end


end
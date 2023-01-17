function [data,meta1,meta2] = fit_channel(T,rbr,channelname,nPeaks,s0,IPoles)


% load('fit_channel.mat')
nSlices = size(T,2);
% mid     = floor(1 + (nSlices-1)/2);
% use     = floor(mid/2);
% id      = (mid - use):(mid + use);

% round 1
%
for i=1:nSlices

    t = subtrbkgrnd(T(:,i),rbr);
    iPoles = IPoles(i,:);

    % side 1
    t1 = ant_to_pos(t,iPoles,1);
    [H,I,bcoor] = find_peak_locations(t1,nPeaks,s0);
    [data1,metadata1] = fit_peaks(t1,channelname,nPeaks,H,I,bcoor,s0);

    % side 2
     t2 = ant_to_pos(t,iPoles,2);
     [H,I,bcoor] = find_peak_locations(t2,nPeaks,s0);
     [data2,metadata2] = fit_peaks(t2,channelname,nPeaks,H,I,bcoor,s0);

    % get single data by combining both sides. Currently, we take the mean.
    % This can change in the future.
    data.sPeaks(i,:)      = (data1.sPeaks + data2.sPeaks)/2;
    data.AntBorder(i,:)   = (data1.AntBorder + data2.AntBorder)/2; 
    data.PosBorder(i,:)   = (data1.PosBorder + data2.PosBorder)/2;
    data.Width(i,:)       = (data1.Width + data2.Width)/2;


    % make meta1 and meta2
    %
    fn = fieldnames(data1);
    for j=1:length(fn)
        meta1.side1.(fn{j}) = data1.(fn{j});
        meta1.side2.(fn{j}) = data2.(fn{j});
    end
    fn = fieldnames(metadata1);
    for j=1:length(fn)
        meta2.side1.(fn{j}) = metadata1.(fn{j});
        meta2.side2.(fn{j}) = metadata2.(fn{j});
    end
 
end





end
function [IM,bg,bgsig] = subtract_background_image(IM,numch)
% clean image by subtracting background
%
% We assume the background intensity level (i.e., the true black level) for 
% each channel is equal to the mode of intensities seen in a slice.
%
bg    = zeros(numch,1);
bgsig = zeros(numch,1);
for i = 1:numch
    slice       = double(IM(:,:,i,1));
    [n,x]       = hist(slice(:),0:4095);
    [A,k]       = max(n(1:end-1));
    bg(i)       = x(k);
    X           = x(max(k-4,1):k+4);
    Y           = n(max(k-4,1):k+4);
    SIG         = (X-bg(i))./sqrt(-2*log(Y/A));
    bgsig(i)    = sqrt(meanDU(SIG'.^2));
    IM(:,:,i,:) = imsubtract(IM(:,:,i,:),bg(i));
end

end


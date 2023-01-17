function [Img, bwlabel] = segment_embryo_image(I)
% This function returns a segmented image of the emrbyo and its label.
% Input:
%           I: Image of any class and a maximum of three dimensions
% Output:
%           Img: Segmented image
%           bwlabel: label, in black and white.

Itmp    = max(I,[],3);              % make image 2D
%se      = strel('disk',2);
%Itmp    = imclose(Itmp,se);         % close image
mg      = mat2gray(Itmp);           % converts image to grayscale
bw      = imbinarize(mg);           % converts inage to bw/binary
bw      = imfill(bw,'holes');       % fills any black regions surrounded by white regions
CC      = bwconncomp(bw);           % finds connected components in the bw image 
CC2     = bwconncomp(imbinarize(mg));

if CC2.NumObjects>1
    rp      = regionprops(CC, 'Area', 'PixelIdxList');
    areas   = zeros(1,CC.NumObjects);
    for i = 1:CC.NumObjects
        areas(i) = rp(i).Area;
    end
    [~, idx] = max(areas);
    Itmp     = zeros(size(bw));
    Itmp(CC.PixelIdxList{idx}) = 1;
    bw       = Itmp;
else
    Itmp     = zeros(size(bw));
    Itmp(CC.PixelIdxList{1})=1;
    bw       = Itmp;
end

% Output
c       = class(I);
Img     = I.*feval(c,bw);           % Final segmented image.
bwlabel = bw;                       % Final label



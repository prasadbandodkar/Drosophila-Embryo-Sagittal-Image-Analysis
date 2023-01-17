function h = make_montage(IM,towrite,location)

if ~exist('towrite','var')
   towrite = false; 
else
    if ~exist('location','var')
        location = 'northwest';
    end
end



% get image dimensions
IM      = squeeze(IM);
[r,c,nIM] = size(IM);
if isodd(nIM)
   n = nIM+1;
else
   n = nIM;
end

% get best rows and columns of the montage
i0          = divisors(n);
i0(i0==1)   = [];
nrc         = length(i0);
if nrc == 1
    rows = 1;
    cols = 1;
else
    i    = round(nrc/2);
    cols = i + 1;
    rows = i;
end

% Initialize final image and text writing vectors
I     = zeros(r*rows,c*cols);
if towrite
    ytext = zeros(nIM,1);
    xtext = ytext;
end
lenr  = 0.05*r;
lenc  = 0.05*c;




% Loop over all images
count    = 1;
istart1  = 1;   istart2  = 1;
iend1    = r;   iend2    = c;  
for i=1:rows
    for j=1:cols
        I(istart1:iend1,istart2:iend2) = mat2gray(IM(:,:,count));    

        if towrite
            switch location
                case 'northwest'
                    ytext(count) = istart1 + lenc;
                    xtext(count) = istart2 + lenr;
                case 'northeast'
                    ytext(count) = istart1 + lenc;
                    xtext(count) = iend2   - lenr;
                case 'southwest'
                    ytext(count) = iend1    - lenr;
                    xtext(count) = istart2  + lenc;            
                case 'southeast'
                    ytext(count) = iend1 - lenc;
                    xtext(count) = iend2 - lenr;
            end
        end

        count = count + 1;
        if count > nIM
           break 
        end
        istart2 = istart2 + c;
        iend2   = iend2 + c;

    end
    istart1 = istart1 + r;
    iend1   = iend1 + r;
    istart2 = 1;
    iend2   = c;
end

dtext = string(1:count-1);
hMontage = figure;
set(hMontage,'visible','off','Position',[1,1,1250,750])  
imshow(I,[])
if towrite
    text(xtext,ytext,dtext,'fontsize',14,'FontWeight','bold','color',[1,1,1])
end
h = gcf;

end
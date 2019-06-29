fileID = fopen('distance.txt', 'wt');
footVid = VideoReader('GOPR0124.mp4');
grayvalue = 50;
sedisk = strel('disk',20);
for iframe = 1:footVid.Duration %* footVid.FrameRate
    grayfoot = rgb2gray(read(footVid,200));
    bwfoot = imextendedmax(grayfoot, grayvalue);
    %imshow(bwfoot);

    noSmallStructures = imopen(bwfoot, sedisk);
    %imshow(noSmallStructures);

    stats = regionprops(noSmallStructures, {'Centroid','Area','BoundingBox'});
   
end
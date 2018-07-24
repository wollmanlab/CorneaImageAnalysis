function [edofimg] = fstackReconst(img,fmap)
%FSTACK merging images of mutiple focal planes into one in-focus image.
%
%  edofimg = fstackReconst(img,fmap) merges img, an img array containing grayscale or
%  color images acquired at mutiple focal distance, into one all-in-focus
%  image.
%  
%  Version 0.1
%  
% Cell array allows more flexible indexing. Some of the logic indexing
% method used here won't work with multi-dimensional array. Logical
% indexing dramatically improves execution speed. 
if ~iscell(img)
    error('Input needs to be cell array of images.')
end


% make brightness of all images equal
avg1 = mean2(img{1});
for ii = 2 : length(img)
    avgcur = mean2(img{ii});
    img{ii} = img{ii} + avg1 - avgcur;
end
    
edofimg = img{1};


% Extract in-focus pixel from every image. 
for ii = 1:length(img)
    index = fmap == ii;
        edofimg(index) = img{ii}(index);
end

% Blend different focal planes
for ii = 1:length(img)-1
    index = fmap > ii & fmap < ii+1;
    edofimg(index) = (fmap(index) - ii).*single(img{ii+1}(index)) + ...
            (ii+1-fmap(index)).*single(img{ii}(index));
end

end


    
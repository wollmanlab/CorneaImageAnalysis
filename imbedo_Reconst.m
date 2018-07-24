function infimg = imbedo_Reconst(Data,fmap1filt2)
infimg = zeros(size(Data,1), size(Data,2));
Nplanes = size(Data,3);

% Extract in-focus pixel from every image. 
for ii = 1:Nplanes
    index = fmap1filt2 == ii;
    tmpimg = Data(:,:,ii);
    infimg(index) = tmpimg(index);
end

% Blend different focal planes
for ii = 1:Nplanes-1
    index = fmap1filt2 > ii & fmap1filt2 < ii+1;
    tmpimg = Data(:,:,ii);
    tmpimg1 = Data(:,:,ii+1);

        infimg(index) = (fmap1filt2(index) - ii).*tmpimg1(index) + ...
            (ii+1-fmap1filt2(index)).*tmpimg(index);
end

end
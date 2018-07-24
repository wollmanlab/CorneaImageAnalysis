function infimg = imbedo_ReconstLMP(Data,fmap1filt2)
infimg = zeros(size(Data,1), size(Data,2));
Nplanes = size(Data,3);


for ii = 1:Nplanes-1
    index = round(fmap1filt2)==ii;% > ii & fmap1filt2 <= ii+1;
    localMaxProj = max(Data(:,:,max(ii-1,1):min(ii+1,Nplanes-1)),[],3);
    infimg(index) = localMaxProj(index);
end

end
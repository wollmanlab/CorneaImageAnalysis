function ind = findPlaneByBrightness(testIMG)
    meanBrightness = squeeze(mean(max(testIMG)));
    ind=find(meanBrightness == max(meanBrightness),1);
    ind(isempty(ind)) = 0;
end
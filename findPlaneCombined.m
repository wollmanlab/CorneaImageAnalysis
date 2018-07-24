function ind = findPlaneCombined(testIMG,lineskip)
    WindowsizeX = size(testIMG,1);
    WindowsizeY = size(testIMG,2);
    testIMG = testIMG(:,(1:floor(WindowsizeY/lineskip))*lineskip,:);
    testIMG = abs(fft(testIMG));
    fftLine2 = squeeze(reshape(testIMG,WindowsizeX*floor(WindowsizeY/lineskip),1,[]));
    fftLine = repmat(abs(fft(ones(WindowsizeX,1))),floor(WindowsizeY/lineskip),1);
    fftLine = repmat(fftLine,1,size(fftLine2,2));
    CorrF = (mean(fftLine.*fftLine2)-mean(fftLine).*mean(fftLine2))./(std(fftLine).*std(fftLine2));
    meanBrightness = squeeze(mean(max(testIMG)));
    CorrF = meanBrightness'./CorrF;
    %findpeaks(-CorrF)
    %pause
    ind=find(CorrF == max(CorrF),1);
    ind(isempty(ind)) = 0;
end
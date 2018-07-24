function [infimg, fmap1filt2] = imbedo(Data,Nbins)

[infimg, fmap1filt2] = imbedo_v2(Data,Nbins);
% Nbins - 2^something
% Data, image stack

% if numel(Nbins)==1
%     Nbins(2) = Nbins(1);
% end;
% WindowsizeX = size(Data,1)/Nbins(1);
% WindowsizeY = size(Data,2)/Nbins(2);
% 
% while(rem(WindowsizeX,1))%input control
%     Nbins(1) = Nbins(1)+1;
%     WindowsizeX = size(Data,1)/Nbins(1);
% end
% 
% while(rem(WindowsizeY,1))%input control
%     Nbins(2) = Nbins(2)+1;
%     WindowsizeY = size(Data,2)/Nbins(2);
% end
% 
% lineskip = 10;
% Nplanes = size(Data,3);
% % make line of ffts for out of focus image
% 
% 
%     fftLine = repmat(abs(fft(ones(WindowsizeX,1))),floor(WindowsizeY/lineskip),1);
% 
% % Find focal planes for all regions
% 
% fmap1 = zeros(Nbins);
% for delX = 0:Nbins(1)-1; %<Nbins
%     for delY =  0:Nbins(2)-1; %<Nbins
%         x0 = WindowsizeX*delX+1;
%         y0 = WindowsizeY*delY+1;
%         testIMG = Data(x0:x0+WindowsizeX-1,y0:y0+WindowsizeY-1,:);
%         %stkshow(testIMG)
%         % Calculate correlation of FFTs
%         CorrF = zeros(Nplanes,1);
%         %CorrF1 = zeros(Nplanes,1);
%         for plane = 1:Nplanes;
%             %fftLine2 = [];
%             fftLine2 = zeros(size(testIMG(:,1,plane),1)*floor(WindowsizeY/lineskip),1);
%             for i=1:floor(WindowsizeY/lineskip);
%                 fftLine2((i-1)*size(testIMG(:,1,plane),1)+1:i*size(testIMG(:,1,plane),1)) =  abs(fft(testIMG(:,i*lineskip,plane)));
%             end;
%             
%             %CorrF(plane) =  corr(fftLine, fftLine2);
%             CorrF(plane) = (mean(fftLine.*fftLine2)-mean(fftLine)*mean(fftLine2))./(std(fftLine)*std(fftLine2));
%         end
%         ind=find(CorrF == min(CorrF),1);
%         ind(isempty(ind)) = 0;
%         fmap1(delX+1, delY+1) = ind;
%          %plot(CorrF)
%          %shg
%          %pause()
%     end
% end
% 
% % dealing with outliers... This is imposing continuity: I find the 5% of pixels that are the least like
% % their nearest neighbors, and replace the focal plane with the average focal plane of
% % the neighbors.
%  %flt = fspecial('average');
%  flt = [1 1 1; 1 0 1; 1 1 1]./8;
% 
%  U = imfilter(fmap1,flt,'replicate');
%  [h,x] = histcounts(reshape((sqrt((fmap1-U).^2)),[],1),100);
%  J = cumsum(h)/sum(h)>.95;
%  ind2switch = find(sqrt((fmap1-U).^2) >= min(x(J)));
%  %
%  for repInd = 1:size(ind2switch)
%      [I,J] = ind2sub([size(fmap1,1) size(fmap1,2)], ind2switch(repInd));
%     if I>1 && I<size(fmap1,1) && J>1 && J<size(fmap1,2)
%      fmap1(I,J) = (fmap1(I+1,J)+fmap1(I-1,J)+fmap1(I,J+1)+fmap1(I,J-1)+...
%                       fmap1(I+1,J+1)+fmap1(I-1,J+1)+fmap1(I+1,J-1)+fmap1(I-1,J-1))/8;                  
%                   %sides
%     elseif I==1 && J>1 && J<size(fmap1,1)
%      fmap1(I,J) = (fmap1(I+1,J)+fmap1(I,J+1)+fmap1(I,J-1)+...
%                       fmap1(I+1,J+1)+fmap1(I+1,J-1))/5;
%     
%     elseif I==size(fmap1,1) && J>1 && J<size(fmap1,2)
%      fmap1(I,J) = (fmap1(I-1,J)+fmap1(I,J+1)+fmap1(I,J-1)+...
%                       fmap1(I-1,J+1)+fmap1(I-1,J-1))/5;        
%     
%     elseif I>1 && I<size(fmap1,1) && J==1
%      fmap1(I,J) = (fmap1(I+1,J)+fmap1(I-1,J)+fmap1(I,J+1)+...
%                       fmap1(I+1,J+1)+fmap1(I-1,J+1))/5;
%     
%     elseif I>1 && I<size(fmap1,1) &&  J==size(fmap1,2)
%      fmap1(I,J) = (fmap1(I+1,J)+fmap1(I-1,J)+fmap1(I,J-1)+...
%                       fmap1(I+1,J-1)+fmap1(I-1,J-1))/5;
%                   %corners
%     elseif I==1 && J==1
%      fmap1(I,J) = (fmap1(I+1,J)+fmap1(I,J+1)+...
%                       fmap1(I+1,J+1))/3;
%     
%     elseif I==size(fmap1,1) && J==1
%      fmap1(I,J) = (fmap1(I-1,J)+fmap1(I,J+1)+...
%                       fmap1(I-1,J+1))/3;        
%     
%     elseif I==1 &&  J==size(fmap1,2)
%      fmap1(I,J) = (fmap1(I+1,J)+fmap1(I,J-1)+...
%                       fmap1(I+1,J-1))/3;
%     
%     elseif I==size(fmap1,1) &&  J==size(fmap1,2)
%      fmap1(I,J) = (fmap1(I-1,J)+fmap1(I,J-1)+...
%                       +fmap1(I-1,J-1))/3;          
%     end;  
%  end
% 
% fmap1 = kron(fmap1, ones(size(Data,1)/Nbins(1), size(Data,2)/Nbins(2)));
% %
% fmap1filt2 = imfilter(fmap1,fspecial('gaussian', [200 200], 30),'replicate');
% fmap1filt2(fmap1filt2 < 1) = 1;
% fmap1filt2(fmap1filt2 > Nplanes) = Nplanes;
% 
% 
% %
% infimg = zeros(size(Data,1), size(Data,2));
% % Extract in-focus pixel from every image. 
% for ii = 1:Nplanes
%     index = fmap1 == ii;
%     tmpimg = Data(:,:,ii);
%     infimg(index) = tmpimg(index);
% end
% %
% % Blend different focal planes
% for ii = 1:Nplanes-1
%     index = fmap1filt2 > ii & fmap1filt2 < ii+1;
%     tmpimg = Data(:,:,ii);
%     tmpimg1 = Data(:,:,ii+1);
% 
%         infimg(index) = (fmap1filt2(index) - ii).*tmpimg1(index) + ...
%             (ii+1-fmap1filt2(index)).*tmpimg(index);
% end

end
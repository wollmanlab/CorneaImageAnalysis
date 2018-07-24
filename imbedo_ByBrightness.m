function [infimg, fmap1filt2] = imbedo_ByBrightness(Data,Nbins)
    
    
    % Nbins - 2^something
    % Data, image stack
    
    if numel(Nbins)==1
        Nbins(2) = Nbins(1);
    end;
    WindowsizeX = size(Data,1)/Nbins(1);
    WindowsizeY = size(Data,2)/Nbins(2);
    
    while(rem(WindowsizeX,1))%input control
        Nbins(1) = Nbins(1)+1;
        WindowsizeX = size(Data,1)/Nbins(1);
    end
    
    while(rem(WindowsizeY,1))%input control
        Nbins(2) = Nbins(2)+1;
        WindowsizeY = size(Data,2)/Nbins(2);
    end
    
    lineskip = 10;
    Nplanes = size(Data,3);
    
    % Find focal planes for all regions
    fmap1 = zeros(Nbins);
    
    %Bless JIT, loop faster than oneliner
    %DataCell = mat2cell(Data,WindowsizeX.*ones(1,Nbins(1)),WindowsizeY.*ones(1,Nbins(2)),Nplanes);

    for delX = 0:Nbins(1)-1; %<Nbins
        for delY =  0:Nbins(2)-1; %<Nbins
            x0 = WindowsizeX*delX+1;
            y0 = WindowsizeY*delY+1;
            testIMG = Data(x0:x0+WindowsizeX-1,y0:y0+WindowsizeY-1,:);
            ind = findPlaneByBrightness(testIMG);
            
            fmap1(delX+1, delY+1) = ind;
        end
    end
    
    flt = [1 1 1; 1 0 1; 1 1 1]./8;
    
    U = imfilter(fmap1,flt,'replicate');
    [h,x] = histcounts(reshape((sqrt((fmap1-U).^2)),[],1),1000);
    J = cumsum(h)/sum(h)>.85;
    ind2switch = find(sqrt((fmap1-U).^2) >= min(x(J)));
    %
    for repInd = 1:size(ind2switch)
        [I,J] = ind2sub([size(fmap1,1) size(fmap1,2)], ind2switch(repInd));
        if I>1 && I<size(fmap1,1) && J>1 && J<size(fmap1,2)
            fmap1(I,J) = (fmap1(I+1,J)+fmap1(I-1,J)+fmap1(I,J+1)+fmap1(I,J-1)+...
                fmap1(I+1,J+1)+fmap1(I-1,J+1)+fmap1(I+1,J-1)+fmap1(I-1,J-1))/8;
            %sides
        elseif I==1 && J>1 && J<size(fmap1,1)
            fmap1(I,J) = (fmap1(I+1,J)+fmap1(I,J+1)+fmap1(I,J-1)+...
                fmap1(I+1,J+1)+fmap1(I+1,J-1))/5;
            
        elseif I==size(fmap1,1) && J>1 && J<size(fmap1,2)
            fmap1(I,J) = (fmap1(I-1,J)+fmap1(I,J+1)+fmap1(I,J-1)+...
                fmap1(I-1,J+1)+fmap1(I-1,J-1))/5;
            
        elseif I>1 && I<size(fmap1,1) && J==1
            fmap1(I,J) = (fmap1(I+1,J)+fmap1(I-1,J)+fmap1(I,J+1)+...
                fmap1(I+1,J+1)+fmap1(I-1,J+1))/5;
            
        elseif I>1 && I<size(fmap1,1) &&  J==size(fmap1,2)
            fmap1(I,J) = (fmap1(I+1,J)+fmap1(I-1,J)+fmap1(I,J-1)+...
                fmap1(I+1,J-1)+fmap1(I-1,J-1))/5;
            %corners
        elseif I==1 && J==1
            fmap1(I,J) = (fmap1(I+1,J)+fmap1(I,J+1)+...
                fmap1(I+1,J+1))/3;
            
        elseif I==size(fmap1,1) && J==1
            fmap1(I,J) = (fmap1(I-1,J)+fmap1(I,J+1)+...
                fmap1(I-1,J+1))/3;
            
        elseif I==1 &&  J==size(fmap1,2)
            fmap1(I,J) = (fmap1(I+1,J)+fmap1(I,J-1)+...
                fmap1(I+1,J-1))/3;
            
        elseif I==size(fmap1,1) &&  J==size(fmap1,2)
            fmap1(I,J) = (fmap1(I-1,J)+fmap1(I,J-1)+...
                +fmap1(I-1,J-1))/3;
        end;
    end
    
    fmap1 = kron(fmap1, ones(size(Data,1)/Nbins(1), size(Data,2)/Nbins(2)));
    %
    %fmap1filt2 = imfilter(fmap1,fspecial('gaussian', [200 200], 30),'replicate');
        fmap1filt2 = imgaussfilt(fmap1,50);
    fmap1filt2(fmap1filt2 < 1) = 1;
    fmap1filt2(fmap1filt2 > Nplanes) = Nplanes;
    
    
    %
    infimg = zeros(size(Data,1), size(Data,2));
    % Extract in-focus pixel from every image.
    for ii = 1:Nplanes
        index = fmap1 == ii;
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

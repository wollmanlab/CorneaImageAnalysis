function CorneaCelllbl = CorneaCellsConstructor(fpath, Well,frame)
    CorneaCelllbl = CorneaCellsLbl;
    CorneaCelllbl.PosName = Well;
    CorneaCelllbl.Frame = frame;
    CorneaCelllbl.pth = fpath;
    
    i=frame;
    MD=Metadata(fpath);
    fpathresults = [fpath filesep 'EDoFandPIVProcessing'];

    PIV_MD = Metadata(fpathresults);
    Data = stkread(MD,'Channel','DeepBlue', 'flatfieldcorrection', false, 'frame', frame, 'Position', Well,'register',false);
    %FocusMap = stkread(PIV_MD,'Type','FocusMap', 'flatfieldcorrection', false, 'frame', frame, 'Position', Well,'register',false);
    %Zplanes = PIV_MD.getSpecificMetadata('Zplanes','Type','FocusMap', 'frame', frame, 'Position', Well);
    %Zplanes = Zplanes{1};

    CorneaCelllbl.ImageDims = size(Data);
    %pixSize = MD.getSpecificMetadata('PixelSize');
    %pixSize =pixSize{1};
    %Zindexes = max(unique(cell2mat(MD.getSpecificMetadata('Zindex'))));
    %dZ = MD.getSpecificMetadata('Z');
    %dZ = dZ{2}-dZ{1};
    zAspect = 11;
    CorneaCelllbl.zAspect = zAspect;
    
    imageSize=CorneaCelllbl.ImageDims;
    %
    Thresh = meanThresh(Data);    
    %
    mask = zeros(imageSize);
    for j=1:size(Data,3)
        mask(:,:,j) = bwconvhull(Data(:,:,j)>Thresh);
    end
    
    % find threshold for cells, and make nice masks
    
    DataFilt = imgaussian3(Data,[5 1]);
    
    candidatePeaks = imregionalmax(DataFilt);
    [xx,yy,zz] = ndgrid(-3:3);
    nhood = sqrt(xx.^2 + yy.^2+zz.^2) <= 3.0;
    ExpandedPeaks = imdilate(candidatePeaks,nhood);
    somePix = DataFilt(ExpandedPeaks);
    %[counts, x] = hist(log(somePix),100);
    %ThreshVal = exp(min(x)+otsuthresh(counts)*(max(x)-min(x)));
    ThreshVal = exp(meanThresh(log(somePix)));
    nuclearSeeds = ExpandedPeaks.*(DataFilt>ThreshVal).*mask;
    nuclearSeedsSignal = Data.*nuclearSeeds;
    
    
    CC = bwconncomp(nuclearSeeds);
    S = regionprops(CC,nuclearSeedsSignal,'MeanIntensity','Centroid','Area');
    % filter CC keep only cells with volume in some range
    Areas = cat(1, S.Area);
    Intensities = cat(1, S.MeanIntensity);
    Centroids = cat(1,S.Centroid);
    %remove specks
    J = Areas>=120;
    Areas = Areas(J);
    Intensities = Intensities(J);
    Centroids = Centroids(J,:);
    S = S(J);
    CC.PixelIdxList = CC.PixelIdxList(J);
    CC.NumObjects = numel(CC.PixelIdxList);
    
    %Correct aspect ratio
    n = size(Centroids,1);
    Centroids = Centroids.*repmat([1 1 zAspect],n,1);
    
    %sort by something reasonable: We will sort the centroids by their
    %topologically corrected hight!
    imageSize = size(Data);
    pad = 0;
    x = -pad+1:imageSize(2)+pad;
    y = -pad+1:imageSize(1)+pad;
    [xx, yy] = meshgrid(x,y);
    
    % 2 rounds to remove outliers
    %indx = datasample(1:n,min(4000,n)); %First, calculate approximate surface with a sample of the points
    %F = scatteredInterpolant(Centroids(indx,1),Centroids(indx,2),Centroids(indx,3),'linear','nearest');
    %zz = F(xx,yy);
    %zz = imgaussfilt(zz,100);
    %fftzz = fft2(zz);%old fashioned low pass filter. Fast!
    %fftzz = fftshift(fftzz);
    %lfp = sqrt((xx-imageSize(2)/2).^2+(yy-imageSize(1)/2).^2)<4;
    %rec = ifft2(lfp.*fftzz);
    %zz = abs(rec);
    F=0;
    [~, zz] = imbedo_ByBrightness(Data.*mask,2^8);
    zz = imgaussfilt(zz,100)*zAspect;


    TopoZ = interp2(xx,yy,zz,Centroids(:,1),Centroids(:,2),'nearest');
    %aha!
    DT = delaunayTriangulation([Centroids(:,1),Centroids(:,2),TopoZ]);
    [~, DistFromManifold] = nearestNeighbor(DT,Centroids);
    DistFromManifold = DistFromManifold.*sign(Centroids(:,3)-TopoZ);
    
    %
    [h,x] = hist(DistFromManifold,1000); %Find ouliers
    cdfdZ = cumsum(h)/sum(h);
    upperdZ = x(find(cdfdZ>0.99,1));
    lowerdZ = x(find(cdfdZ>0.02,1)-1);
    indx = logical((DistFromManifold>lowerdZ).*(DistFromManifold<upperdZ));%indices to keeo for final calculation
    
%     F = scatteredInterpolant(Centroids(indx,1),Centroids(indx,2),Centroids(indx,3),'linear','nearest');
%     zz = F(xx,yy);
%     zz = imgaussfilt(zz,100);
% 
%     TopoZ = interp2(xx,yy,zz,Centroids(:,1),Centroids(:,2),'nearest');
%     %aha!
%     DT = delaunayTriangulation([Centroids(:,1),Centroids(:,2),TopoZ]);
%     [~, DistFromManifold] = nearestNeighbor(DT,Centroids);
%     DistFromManifold = DistFromManifold.*sign(Centroids(:,3)-TopoZ);
%     
%     %
%     [h,x] = hist(DistFromManifold,1000); %Find ouliers
%     cdfdZ = cumsum(h)/sum(h);
%     upperdZ = x(find(cdfdZ>0.99,1));
%     lowerdZ = x(find(cdfdZ>0.03,1)-1);
%     indx = logical((DistFromManifold>lowerdZ).*(DistFromManifold<upperdZ));%indices to keeo for final calculation
    
    Centroids = Centroids(indx,:);
    TopoZ = TopoZ(indx);
    DistFromManifold = DistFromManifold(indx);
    Intensities = Intensities(indx);
    Areas=Areas(indx);
    S = S(indx);
    CC.PixelIdxList = CC.PixelIdxList(indx);
    CC.NumObjects = numel(CC.PixelIdxList);
    n = size(Centroids,1);
    
    
    [~,J] = sort(DistFromManifold); %Sort by distance from smoothed topological map
    TopoZ = TopoZ(J);
    Centroids = Centroids(J,:);
    Intensities = Intensities(J);
    CC.PixelIdxList = CC.PixelIdxList(J);
    
    
    %Apply drift correction to centroids!
    Tforms = MD.getSpecificMetadata('driftTform','Position',Well, 'frame', i);
    dY = Tforms{1}(7);
    dX = Tforms{1}(8);
    Centroids = Centroids+repmat([dY, dX, 0],n,1);
    CorneaCelllbl.Centroids = Centroids;
    CorneaCelllbl.Intensities = Intensities;
    CorneaCelllbl.num = n;
    %CorneaCelllbl.DistFromManifold = DistFromManifold;
    CorneaCelllbl.TopoZ = TopoZ;
    %CorneaCelllbl.GridZ = zz;
    CorneaCelllbl.CC = CC;
    i
    
    
    
end

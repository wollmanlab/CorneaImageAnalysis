function [Wound, Masks] = WoundTracker(Red)
R1 = Red(:,:,1); %look at first Red layer
[h,x] = hist(R1(:),100); %make histogram
mask1 = R1>x(find(cumsum(h)./sum(h)>0.95,1,'first'));%mask top 5% red pixels
figure('name','draw around wound')
imagesc(mask1)
shg
pol1 = impoly;
mask1 = logical(mask1.*pol1.createMask);
imagesc(mask1)
shg

% NumObj = 5;
% filtSize = 0;
% while NumObj>1 && filtSize<=20;
% filtSize = filtSize+1
% se = strel('disk',1);
% closeBW = imclose(gpuArray(mask1),se);
% se = strel('disk',filtSize);
% closeBW = imerode(gpuArray(closeBW),se);
% se = strel('disk',filtSize);
% closeBW = imdilate(gpuArray(closeBW),se);
% CC = bwconncomp(gather(closeBW));
% NumObj =  CC.NumObjects;
% imagesc(closeBW);shg;
% pause
% end
se = strel('disk',5);
closeBW = imerode(mask1,se);
filtSize = 1;
InitWoundMask = bwconvhull(closeBW);
se = strel('disk',filtSize);
InitWoundMask = imerode(InitWoundMask,se);

 imagesc(R1.*InitWoundMask,[0.05 0.08])
 InitWound = R1.*InitWoundMask;
 
 CC = bwconncomp(InitWoundMask);
 S0 = regionprops(CC,'BoundingBox','Area','Centroid');

 WoundCenterXY0 = ceil([S0.BoundingBox(1)+S0.BoundingBox(3)/2 S0.BoundingBox(2)+S0.BoundingBox(4)/2]);
 
axis equal;
ROI = ceil(S0.BoundingBox+[-100 -100 200 200])

movROI = R1(ROI(2)+1:ROI(2)+ROI(4),ROI(1)+1:ROI(1)+ROI(3));figure;
imagesc(movROI);shg;

%% Within a ROI, we look for the best registration map. we apply this map to the mask, and move the ROI to the center of the mask
NextROI = ROI;
NextWoundMask = InitWoundMask;
E = edge(NextWoundMask,'canny');
E = imdilate(E,strel('disk',30));
E = bwmorph(E,'thin',inf);

NextWoundMask0 = bwconvhull(E);
NextWoundMask = NextWoundMask0;

WoundCenterXY = WoundCenterXY0;

% In next slice, look at ROI.
close all
Wound=[];

Wound(1).Mask(:,:) =  NextWoundMask;
Masks(:,:,1) = NextWoundMask;
Wound(1).XY = WoundCenterXY;
Wound(1).Centroid = S0.Centroid;

for i=1:size(Red,3)-1;
sprintf('%d out of %d',i,size(Red,3)-1)
R1 = Red(:,:,i);
R2 = Red(:,:,i+1); %look at first Red layer
NWMROI = NextWoundMask(NextROI(2)+1:NextROI(2)+NextROI(4),NextROI(1)+1:NextROI(1)+NextROI(3));
movROI = (R1(NextROI(2)+1:NextROI(2)+NextROI(4),NextROI(1)+1:NextROI(1)+NextROI(3)));
fixROI = (R2(NextROI(2)+1:NextROI(2)+NextROI(4),NextROI(1)+1:NextROI(1)+NextROI(3)));
[Tform,~]= imregdemons(movROI,fixROI,'AccumulatedFieldSmoothing',2,'PyramidLevels',6);
NWMROI = imwarp(NWMROI,(Tform));

NextWoundMask(NextROI(2)+1:NextROI(2)+NextROI(4),NextROI(1)+1:NextROI(1)+NextROI(3)) = NWMROI(1:NextROI(4),1:NextROI(3));


E = edge(NextWoundMask,'canny');
E = imdilate(E,strel('disk',30));
E = bwmorph(E,'thin',inf);

NextWoundMask = bwconvhull(E);
CC = bwconncomp(NextWoundMask);
S = regionprops(CC,'BoundingBox','Centroid','Area');
S = S([S.Area]>0.5*S0.Area);


while S.Area>1.1*S0.Area; %keeping size in control
    NextWoundMask = imerode(NextWoundMask,strel('disk',1));
    CC = bwconncomp(NextWoundMask);
    S = regionprops(CC,'BoundingBox','Area','Centroid');
end


WoundCenterXY = ceil([S.BoundingBox(1)+S.BoundingBox(3)/2 S.BoundingBox(2)+S.BoundingBox(4)/2]);
NextROI = ceil(S.BoundingBox+[-100 -100 200 200])
if NextROI(1)<0
    NextROI(3) = NextROI(3)+NextROI(1);
    NextROI(1) = 0;
elseif NextROI(1)+NextROI(3)>size(Red,1)
    NextROI(3) = size(Red,1)-NextROI(1);
end

if NextROI(2)<0
    NextROI(4) = NextROI(4)+NextROI(2);
    NextROI(2) = 0;
elseif NextROI(2)+NextROI(4)>size(Red,2)
    NextROI(4) = size(Red,2)-NextROI(2);
end

Wound(i+1).Mask(:,:) =  NextWoundMask;
Masks(:,:,i+1) = NextWoundMask;
Wound(i+1).XY = WoundCenterXY;
Wound(i+1).Centroid = S.Centroid;
end

end
%% Load results
%% Single cell analysis for incuscope cornea data
BaseStr = regexprep([char(ispc.*'Z:\Images2018\') char(isunix.*'/bigstore/Images2018/')],char(0),'');
Usr = 'Jen';
Project = 'CorneaCCM';
Dataset = 'WoundAgarTitr_2018Mar16_2018May01';
acquisition = 2;
%% Get MD of SAR data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname filesep 'SAR'];
MD=Metadata(fpath);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
Zindexes = max(unique(cell2mat(MD.getSpecificMetadata('Zindex'))));
zAspect = 11;
%%
R = MultiPositionSingleCellWoundResults(fpath)
%%
%filename = '/bigstore/GeneralStorage/Alon/figsForLM/HeatMapsRadialVelocity.gif';
fpath = '/bigstore/GeneralStorage/Alon/figsForLM/';
outputVideo = VideoWriter(sprintf('%sHeatMapsSpeed%s',fpath,'.avi'));
outputVideo.FrameRate = 7;
open(outputVideo)

n=1
for i=1:220
R.HeatMapData('Speed',R.PosNames{1},i,'clims',[0,15])
axis tight
set(gcf,'color','k')
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
title('Speed','color','w')
ch = get(gcf,'Children');
set(ch,'color','w')
PicComment_v1('str', [num2str((i-1)*.5) 'h'],'fontsize', 10,'color','g')
drawnow
ColorbarLabel = get(ch(1),'Label');
ColorbarLabel.String = '\mum / hour';



frame = getframe(gcf);
im = frame2im(frame);
%[imind,cm] = rgb2ind(im,256);
% if n == 1
%     imwrite(imind,cm,filename,'gif', 'Loopcount',1,'DelayTime',.20);
% else
%     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.20);
% end
% n=n+1;
   writeVideo(outputVideo,im);

end
close(outputVideo)
%%
j=1
CorneaCells = R.getCorneaCellsLbl(R.PosNames{j});
CorneaCells{7}.scatter3('epi');
set(gcf,'color','k')
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gca,'zcolor','w')
      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'ScatterEpi']);

%% Example of peeling layers
h = figure('Position',[100 100 1000 500],'color','k')
axis tight manual
filename = '/bigstore/GeneralStorage/Alon/figsForLM/AnimatedLayerPeeling.gif';

j=1
CorneaCells = R.getCorneaCellsLbl(R.PosNames{j});

n=1;
i=1

for j=0:20
CorneaCells{i}.scatter3('epi',find(CorneaCells{i}.epiScore>0.5),25*j)
hold on;
CorneaCells{i}.scatter3('epi',find(CorneaCells{i}.epiScore<0.5))
hold off;
set(gca,'Position',[0.3, 0.1, 0.4, 0.8])
set(gcf,'color','k')
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gca,'zcolor','w')
drawnow
shg
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if n == 1
    imwrite(imind,cm,filename,'gif', 'Loopcount',1,'DelayTime',.20);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.20);
end
n=n+1;
end

%
for k=1:10
CorneaCells{i}.scatter3('epi',find(CorneaCells{i}.epiScore>0.5),25*j)
hold on;
CorneaCells{i}.scatter3('epi',find(CorneaCells{i}.epiScore<0.5))
hold off;
set(gca,'Position',[0.3, 0.1, 0.4, 0.8])
set(gcf,'color','k')
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gca,'zcolor','w')
shg
 drawnow
 pause
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if n == 1
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',.20);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.20);
end
n=n+1;
end



for i=51:200
    clf
    ax1 = axes('Position',[0.1, 0.1, 0.4, 0.8])
CorneaCells{i}.scatter('epi',find(CorneaCells{i}.epiScore>0.5),25*j)
set(ax1,'xdir','reverse','ydir','reverse')
set(gcf,'color','k')
set(ax1,'color','k')
set(ax1,'xcolor','w')
set(ax1,'ycolor','w')
title('Upper epithelium','color','w')
    ax2 = axes('Position',[0.6, 0.1, 0.4, 0.8])
CorneaCells{i}.scatter('epi',find(CorneaCells{i}.epiScore<0.5))
set(ax2,'xdir','reverse','ydir','reverse')
set(ax2,'color','k')
set(ax2,'xcolor','w')
set(ax2,'ycolor','w')
title('Lower epithelium','color','w')

drawnow
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if n == 1
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',.20);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.20);
end
n=n+1;
end








%%
a = R.allTrackMatrix(R.PosNames{1});
imagesc(a)
set(gca,'xtick',[1:50:251],'xticklabel',[0:50:250]/2)
colormap(magma)
xlabel('Time(h)')
ylabel('Track#')
set(gcf,'color','k')
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'AllTrackMat']);



%%
plotCompTrack(R.Tracks{1}(2))
set(gcf,'color','k')
ch = get(gcf,'Children');
set(ch,'color','k')
set(ch,'xcolor','w')
set(ch,'ycolor','w')
      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'SingleTrackExample1']);
%%
plotTracks2D(R.Tracks{1}(2),[],1)
set(gcf,'color','k')
ch = get(gcf,'Children');
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gca,'xlim',[2000    2433],'ylim',[1497    2000])
ch(1).delete;
ch2 = get(ch(2),'Children')
ch2(1).CData=[1 1 1]
ch2(2).CData=[1 1 1]
      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'SingleTrackPathExample1']);
%%
plotCompTrack(R.Tracks{1}(100))
set(gcf,'color','k')
ch = get(gcf,'Children');
set(ch,'color','k')
set(ch,'xcolor','w')
set(ch,'ycolor','w')
      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'SingleTrackExample2']);

%%
plotTracks2D(R.Tracks{1}(100),[],1)
set(gcf,'color','k')
ch = get(gcf,'Children');
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gca,'xlim',[2155    2330],'ylim',[ 1405    1686])
ch(1).delete;
ch2 = get(ch(2),'Children')
ch2(1).CData=[1 1 1]
ch2(2).CData=[1 1 1]
      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'SingleTrackPathExample2']);

%%
J = logical(([R.Tracks{1}.epiRank]<0.52).*([R.Tracks{1}.epiRank]>0.50));
plotTracks2D(R.Tracks{1}(J),[],1);
set(gcf,'color','k')
ch = get(gcf,'Children');
ch(1).delete;
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
ch2 = get(ch(2),'Children')
ch2(1).CData=[1 1 1]
ch2(2).CData=[1 1 1]
      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'ManyTracksExample']);

%%
figure
J = logical(([R.Tracks{1}.epiRank]<0.05).*([R.Tracks{1}.epiRank]>0));
J = find(J);
for i=1:numel(J)
    
    T = R.Tracks{1}(J(i)).T(1:end);
    D = R.Tracks{1}(J(i)).DistFromWound;
    nanInds = find(isnan(D))
    if any(nanInds)%fill gaps with linear interpolation and dashed lines
        G = nan(size(D));
        skipInds = cumsum([1; diff(find(isnan(D)))>1]);
        
        for ind=1:numel(unique(skipInds))
            jx = skipInds==ind;s
            indsToFill = nanInds(jx);
            m = (D(indsToFill(end)+1)-D(indsToFill(1)-1))/(1+sum(jx));
            G(indsToFill(1)-1:indsToFill(end)+1) = D(indsToFill(1)-1)+m*(0:sum(jx)+1);
        end
        h = plot(T,D,'-',T,G,'-.');
        set(h(2),'color',get(h(1),'color'));
    else
        h = plot(T,D,'-');
    end
    hold all
end
shg
set(gcf,'color','k')
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
title('Top 5%','color','w')
xlabel('frame','color','w')
ylabel('Distance From Wound','color','w')
      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'DistTracksTop']);
%%
J = logical(([R.Tracks{1}.epiRank]<0.05).*([R.Tracks{1}.epiRank]>0.0));
JDist = arrayfun(@(x) x.DistFromWound(1), R.Tracks{1})<250;
JlowClose = logical(JDist.*J');


J = find(JlowClose);
for i=1:numel(J)
plot(R.Tracks{1}(J(i)).T(1:end), R.Tracks{1}(J(i)).DistFromWound);
hold all
end
shg
set(gcf,'color','k')
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
title('Top 5%','color','w')
xlabel('frame','color','w')
ylabel('Distance From Wound','color','w')
      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
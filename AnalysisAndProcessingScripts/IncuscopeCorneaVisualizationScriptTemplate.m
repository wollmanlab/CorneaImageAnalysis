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
%% Load results object
R = MultiPositionSingleCellWoundResults(fpath)
%% Analysis of PIV Data
%% Heatmaps of data over time
%filename = '/bigstore/GeneralStorage/Alon/figsForLM/HeatMapsRadialVelocity.gif';
fpath = '/bigstore/GeneralStorage/Alon/FiguresIncuWound090618/';
outputVideo = VideoWriter(sprintf('%sHeatMapsSpeed%s',fpath,'.avi'));
outputVideo.FrameRate = 7;
open(outputVideo)

n=1
for i=1:220
    %Change the data type for other maps
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


%% Single Cell stuff
%% Example of single epithelium

j=1e
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

for j=0:20 %separate layers
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
for k=1:10 %Rotate manually
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



for i=1:200 %evolve in time independently
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


%% Density Maps
j = 4;
CorneaCells = R.getCorneaCellsLbl(R.PosNames{j});
outputVideo = VideoWriter(sprintf('%sDensityMapsBottom%s%s',fpath,R.PosNames{j},'.avi'));
outputVideo.FrameRate = 7;
open(outputVideo)
for i=1:numel(CorneaCells)
    ccellLbl = CorneaCells{i};
    xx = ccellLbl.Centroids(:,1);
    yy = ccellLbl.Centroids(:,2);
    J = find(ccellLbl.epiScore>=0.65);
    
    [h ,Density,Bins, Color] = PlotColorCoded_v2(xx,yy,J, J, range(xx)/40, range(yy)/40 );shg
    set(gcf,'color','k')
    set(gca,'Position',[0 0 1 1],'xlim',[0 3000],'ylim',[0 2100],'color','k','CLim', [0,2])
    drawnow;
    
    
    frame = getframe(gcf);
    im = frame2im(frame);
    writeVideo(outputVideo,im);
    
end
close(outputVideo)



%% Speed maps
figure()
    set(gcf,'color','k')
j=3
Tracks = R.getTracks(R.PosNames{j});
JinEpi = [Tracks.epiRank]>0.65;
Tracks = Tracks(JinEpi);

outputVideo = VideoWriter(sprintf('%sRadVMapsBottom%s%s',fpath,R.PosNames{j},'.avi'));
outputVideo.FrameRate = 7;
open(outputVideo)
for i=1:numel(R.Frames)-1; %timepoint
    trackNumsThatPassThroughi = find(arrayfun(@(x) any(find(x.T==i).*find(x.T==i+1)), Tracks));
    Speeds = zeros(1,numel(trackNumsThatPassThroughi));
    xx = zeros(1,numel(trackNumsThatPassThroughi));
    yy = zeros(1,numel(trackNumsThatPassThroughi));
    
    for ind1=1:numel(trackNumsThatPassThroughi)
        trackNum = trackNumsThatPassThroughi(ind1);
        track1 = Tracks(trackNum);
        
        indWhenPassed = find(track1.T==i);
        Vel = (track1.RadialVelocity);
        Speeds(ind1) =  Vel(indWhenPassed);
        
        xx(ind1) =  track1.tracksCoordAmpCG(1+8*indWhenPassed);
        yy(ind1) =  track1.tracksCoordAmpCG(2+8*indWhenPassed);
        
    end
    %[~,Jplot] = sort(Speeds);
    h = scatter(xx, yy,[],Speeds,'o','filled');%AOY

    set(gca,'Position',[0 0 1 1],'xlim',[0 3000],'ylim',[0 2100],'color','k','clim',[-10,10])
    axis equal
    colormap(makeColorMap([0.8 0.8 0], [0 0 0], [.8,0,0.6]))
    set(h,'MarkerFaceAlpha', 0.7)
    
     drawnow;
   
    frame = getframe(gcf);
    im = frame2im(frame);
    writeVideo(outputVideo,im);
end
close(outputVideo)


%% Tracks
%% Matrix of all of the tracks
a = R.allTrackMatrix(R.PosNames{1});
imagesc(a)
set(gca,'xtick',[1:50:251],'xticklabel',[0:50:250]/2)
colormap(magma)
xlabel('Time(h)')
ylabel('Track #')
set(gcf,'color','k')
set(gca,'color','k')
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'AllTrackMat']);



%% Example of a single tracked cell

plotCompTrack(R.Tracks{1}(2))
set(gcf,'color','k')
ch = get(gcf,'Children');
set(ch,'color','k')
set(ch,'xcolor','w')
set(ch,'ycolor','w')
%      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
%print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'SingleTrackExample1']);
%% track
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
%      set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
%print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'SingleTrackPathExample1']);

%% Spaghetti bowl
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

%% plot track properties over time
figure
J = logical(([R.Tracks{1}.epiRank]<0.1).*([R.Tracks{1}.epiRank]>0));
J = find(J);
for i=1:numel(J)
    
    T = R.Tracks{1}(J(i)).T(1:end);
    %Change "DistFromWound to whatever else you want to plot
    D = R.Tracks{1}(J(i)).DistFromWound;
    nanInds = find(isnan(D))
    if any(nanInds)%fill gaps with linear interpolation and dashed lines
        G = nan(size(D));
        skipInds = cumsum([1; diff(find(isnan(D)))>1]);
        
        for ind=1:numel(unique(skipInds))
            jx = skipInds==ind;
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
%%
%set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
%print(gcf,'-dpng','-r300',['/bigstore/GeneralStorage/Alon/figsForLM/' 'DistTracksTop']);
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









%% Speeds

j=4
Tracks = R.getTracks(R.PosNames{j});

JinEpi = [Tracks.epiRank]>0.65;
Tracks = Tracks(JinEpi)
%xCorrs = {};
%xCorrsErr = {};
%Ts = {};

T = 1:numel(R.TimeVecs{1})-1;
%T = [-fliplr(T(2:end)),T];
Speeds = nan(numel(Tracks), numel(T));
Dists = nan(numel(Tracks), numel(T));

for ind = 1:numel(Tracks)
    ind

    a = (Tracks(ind).RadialVelocity');
    b = (Tracks(ind).DistFromWound(1:end-1)');
    
    T1 = Tracks(ind).T;
    T1 = T1(1:end-1);

    
    J = ismember(T,T1);
    Speeds(ind,J) = a;
    Dists(ind,J) = b;
    
end





%% V/D
figure('color','w')
Twin = 3;
Start = 1
for i=0:50
    
    J = find(~isnan(nanmean(Speeds(:,Start-1 + Twin*i+[1:Twin]),2)))
    [Xb, Yb, stdXb, stdYb, steXb, steYb] = BinData_v1(nanmean(Dists(J,Start-1 + Twin*i+[1:Twin]),2), nanmean(Speeds(J,Start-1 +Twin*i+[1:Twin]),2), 10);
    ploterr(Xb*PixelSize,Yb,steXb,steYb);
    hold all
    shg
    xlabel('Distance From Wound')
    ylabel('Radial Velocity')
    %pause
end
set(gca,'xtick',[30:60:600]/PixelSize, 'xticklabel',30:60:600)
set(gca,'ylim',[-5,5],'xlim',[0,600])

Ch = get(gca,'Children');
tzeva = flipud(plasma(size(Ch,1)/3));
for i=0:size(Ch,1)/3-1
    for jj=1:3
        Ch(3*i+jj).Color = tzeva(i+1,:);
        Ch(3*i+jj).Color(4) = 0.4;
        Ch(3*i+jj).LineWidth = 1.5
    end
end

cb = colorbar(gca)
cb.Direction = 'reverse'
colormap(flipud(plasma(size(Ch,1))))
cb.Ticks = 0:(1/11):1
cb.TickLabels = 110:-10:0
xlabel(cb,'Time (h)')


set(gcf, 'PaperPositionMode','auto')

print(gcf,'-dpng','-r300',sprintf('%sVelocityVsDistanceBottom%s',fpath,R.PosNames{j}));
%% V/T at D
figure

for D = 50:25:950;
    
    inds = find((Dists<D+50).*(Dists>D));
    [~,iy] = ind2sub(size(Dists),inds);
    S1 = Speeds(inds);
    
    [Xb, Yb, stdXb, stdYb, steXb, steYb] = BinData_v1(iy, S1, 220);
    ploterr(Xb,Yb,steXb,steYb,'-');
    hold all;
end

set(gca,'ylim',[-10,10],'xlim',[0,220],'xtick',[0:20:220], 'xticklabel',[0:20:220]/2)
shg
xlabel('Time(h)')
ylabel('Radial Velocity (\mum/h)')


Ch = get(gca,'Children');
tzeva = flipud(plasma(size(Ch,1)/3));
for i=0:size(Ch,1)/3-1
    for jj=1:3
        Ch(3*i+jj).Color = tzeva(i+1,:);
        Ch(3*i+jj).Color(4) = 0.4;
        Ch(3*i+jj).LineWidth = 1.5
    end
end
cb = colorbar(gca)
cb.Direction = 'reverse'
colormap(flipud(plasma(size(Ch,1))))
cb.Ticks = 0:1/9:1
cb.TickLabels = [570:-60:30]
xlabel(cb,'Distance From Wound (\mum)')


set(gcf, 'PaperPositionMode','auto')

print(gcf,'-dpng','-r300',sprintf('%sVelocityVsTimeBottom%s',fpath,R.PosNames{j}));

%% V/D at T ribbons
figure
Ybs = [];
Xbs = [];
for D = 50:25:900;
    
    inds = find((Dists<D+50).*(Dists>D));
    [~,iy] = ind2sub(size(Dists),inds);
    S1 = Speeds(inds);
    
    [Xb, Yb, stdXb, stdYb, steXb, steYb] = BinData_v1(iy, S1, 220);
    Xbs = [Xbs ; Xb];
    Ybs = [Ybs ; Yb];
    
    hold all;
end

h = ribbon(Xbs', Ybs')
for i=1:numel(h)
    h(i).FaceAlpha=0.6;
end
ylabel('Time (h)')
xlabel('Distance from wound (\mum)')
zlabel('Radial Velocity (\mum/h)')
set(gca,'ytick',[0:40:220], 'yticklabel',[0:40:220]/2,'xtick',25\[30:60:600]/PixelSize, 'xticklabel',30:60:600,'CameraPosition', [-216.5420 -1.3176e+03 40.1221],'zlim',[-5, 5])
colormap((plasma))

set(gcf, 'PaperPositionMode','auto')

print(gcf,'-dpng','-r300',sprintf('%sVelocityVsTimeRibbonsBottom%s',fpath,R.PosNames{j}));

%% V/D at T heatmap
figure
imagesc(Ybs)
colormap(colormap(makeColorMap([0.8 0.8 0], [0 0 0], [.8,0,0.8])))
set(gca,'clim',[-5,5])
cb = colorbar
xlabel('T (h)')
ylabel('Distance from wound (\mum)')
set(gca,'xtick',[0:20:220], 'xticklabel',[0:20:220]/2,'ytick',25\[30:60:600]/PixelSize, 'yticklabel',30:60:600)
xlabel(cb,'Radial Velocity')

set(gcf, 'PaperPositionMode','auto')

print(gcf,'-dpng','-r300',sprintf('%sVelocityVsTimeHeatMapBottom%s',fpath,R.PosNames{j}));

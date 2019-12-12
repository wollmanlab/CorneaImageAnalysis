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



%% Create a results object
R = MultiPositionSingleCellWoundResults();
R.pth = fpath;
R.PosNames=unique(MD.getSpecificMetadata('Position'));

frames = unique(cell2mat(MD.getSpecificMetadata('frame')));

R.Frames = frames;
R.analysisScript=fullfile([fpath filesep 'AnalysisScriptTemplate.m']);%Change this to the right file
R.reportPth = [BaseStr 'Reports' filesep 'Alon' filesep Project filesep Dataset];
Wells = R.PosNames;


%% Add PIV to results, mask wounds, etc.
fpathresults = [MD.pth filesep 'EDoFandPIVProcessing'];

for WellNum=1:numel(Wells)
    pos = Wells{WellNum};
    R.setPIVLbl(PIVConstructor(fpathresults, pos),pos)
end

%%
fpathresults = [fpath filesep 'EDoFandPIVProcessing']
PIV_MD = Metadata(fpathresults);
Wells = unique(PIV_MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(PIV_MD.getSpecificMetadata('frame')));



%% Make basic measurements
dt = 24*(R.TimeVecs{1}(2)-R.TimeVecs{1}(1));
PixelSize = MD.unique('PixelSize');
for j=1:numel(R.PosNames);
    PIV = R.getPIVLbl(R.PosNames{j});
    DistanceFromWound = zeros(size(PIV.X,1),numel(PIV.flist));
    RadialVelocity = zeros(size(PIV.X,1),numel(PIV.flist));
    TangentVelocity = zeros(size(PIV.X,1),numel(PIV.flist));
    Speed = zeros(size(PIV.X,1),numel(PIV.flist));
    
    Cent = PIV.WoundLbl.SmoothCentroid(1:numel(PIV.flist)+1);
    %Drift transforms
    %Tforms = cell2mat(PIV_MD.getSpecificMetadata('driftTform','Type','EDOF', 'Channel','DeepBlue', 'Position', Wells{j}));
    %Tforms = Tforms(:,7:8);
    for i=1:numel(PIV.flist)
        sprintf('frame %d out of %d on position %d out of %d',i,numel(PIV.flist),j,numel(R.PosNames))
        
        x = PIV.X -Cent(i,1);
        y = PIV.Y -Cent(i,2);

        CosAlpha = x./(sqrt(x.^2+y.^2));
        CosAlpha(isnan(CosAlpha))=0;
        SinAlpha = y./(sqrt(x.^2+y.^2));
        SinAlpha(isnan(SinAlpha))=0;
        
        DistanceFromWound(:,i) = PixelSize*sqrt(x.^2+y.^2);
        UMasked = PixelSize*PIV.UMasked(i)/dt;%drift correction on u and v
        VMasked = PixelSize*PIV.VMasked(i)/dt;%all in micron/h
        
        woundU = 0;%mean(UMasked);
        woundV = 0;%mean(VMasked);
        
        RadialVelocity(:,i) = (CosAlpha.*(UMasked-woundU)+SinAlpha.*(VMasked-woundV)).*~PIV.indInWound(i);
        TangentVelocity(:,i) = (SinAlpha.*(UMasked-woundU)-CosAlpha.*(VMasked-woundV)).*~PIV.indInWound(i);
        %UMasked = UMasked+Tforms(i+1,1)-Tforms(i,1);%drift correction on u and v
        %VMasked = VMasked+Tforms(i+1,2)-Tforms(i,2);
        Speed(:,i) = sqrt(UMasked.^2+VMasked.^2).*~PIV.indInWound(i);
    end
    R.setData('DistanceFromWound',DistanceFromWound, R.PosNames{j})
    R.setData('RadialVelocity',RadialVelocity, R.PosNames{j})
    R.setData('TangentVelocity',TangentVelocity, R.PosNames{j})
    R.setData('Speed',Speed, R.PosNames{j})
end
%%
R.saveResults

%% Add CorneaCellsLbl to results, to single cell segmentation
for WellNum=1:numel(Wells)
    pos = Wells{WellNum};
    CorneaCelllbl = cell(numel(frames),1);
    parfor i=R.Frames'
        CorneaCelllbl{i} = CorneaCellsConstructor(fpath, pos,i)
    end
    R.setCorneaCellsLbl(CorneaCelllbl,pos)
end

%% Find epithelium and generate epithelium related features
for i=1%:numel(R.PosNames)
    R.setEpitheliumFitParams(R.PosNames{i})
end

%% refine 
for i=1:numel(R.PosNames)
    i
    R.refineEpitheliumSurface(R.PosNames{i})
end
%% repeat find
for i=1:numel(R.PosNames)
    R.setEpitheliumFitParams(R.PosNames{i})
end
%% Link adjecent frames
for WellNum=1:numel(R.PosNames)  
    R.Link(R.PosNames(WellNum))
end
%% Close gaps
for WellNum=1:numel(Wells)
    R.closeGaps(R.PosNames(WellNum))
end
%%
fpathresults = [fpath filesep 'EDoFandPIVProcessing']
PIV_MD = Metadata(fpathresults);

%% Calculate more features about single cell tracks
PixelSize = MD.unique('PixelSize');
for WellNum=1:numel(Wells)
    Tracks = R.getTracks(R.PosNames{WellNum}); 
    
    %Time
    T = arrayfun(@(x) (x.seqOfEvents(1,1):x.seqOfEvents(2,1)), Tracks,'UniformOutput',false);
    [Tracks.('T')] = T{:};
    %RMS displacement to cells
    RMSDisp = arrayfun(@(x) PixelSize*createDistanceMatrix([x.tracksCoordAmpCG(1)',x.tracksCoordAmpCG(2)',x.tracksCoordAmpCG(3)'],[x.tracksCoordAmpCG(1:8:end)',x.tracksCoordAmpCG(2:8:end)',x.tracksCoordAmpCG(3:8:end)']), Tracks,'UniformOutput',false);
    [Tracks.('RMSDisp')] = RMSDisp{:};
    
    %epiScore within each track, i.e., position within the epithelium
    CorneaCells = R.getCorneaCellsLbl(R.PosNames{WellNum});
    Densities = cellfun(@(x) x.Density(50), CorneaCells,'UniformOutput', false);

    for i=1:numel(Tracks)
        i
        epiScoreTrack{i} = zeros(1,numel(Tracks(i).T));
        DensityTrack{i} = zeros(1,numel(Tracks(i).T));
        
        for j=1:numel(Tracks(i).T)
            if Tracks(i).tracksFeatIndxCG(j)
                epiScoreTrack{i}(j) = CorneaCells{Tracks(i).T(j)}.epiScore(Tracks(i).tracksFeatIndxCG(j));
                DensityTrack{i}(j) = Densities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
            else
                epiScoreTrack{i}(j) = NaN;
                DensityTrack{i}(j) = NaN;
            end
        end
    end
    [Tracks.('epiScoreTrack')] = epiScoreTrack{:};
    [Tracks.('DensityTrack')] = DensityTrack{:};
    
    
    meanEpiScore = arrayfun(@(x) nanmean(x.epiScoreTrack),Tracks,'UniformOutput',false);
    [Tracks.('meanEpiScore')] = meanEpiScore{:};
    
    
    [h,x] = hist([meanEpiScore{:}]',250);
    dataCDF = cumsum(h)/sum(h);
    epiRank = interp1(x,dataCDF,[meanEpiScore{:}]');
    epiRank(isnan([meanEpiScore{:}]'))=NaN;
    for i=1:numel(Tracks)
        Tracks(i).epiRank = epiRank(i);
    end
    
    % Distance from wound
    PIVlbl = R.getPIVLbl(R.PosNames{WellNum});
    WoundCentroidXY = PIVlbl.WoundLbl.SmoothCentroid;
    Disps = arrayfun(@(x) PixelSize*([x.tracksCoordAmpCG(1:8:end)',x.tracksCoordAmpCG(2:8:end)']-WoundCentroidXY(x.T,:)),Tracks,'UniformOutput',false);
    [Tracks.('DisplacementFromWound')] = Disps{:};
        Disps = cellfun(@(x) [Smoothing(x(:,1)),Smoothing(x(:,2))],Disps,'UniformOutput',false);

    Dists = cellfun(@(x) sqrt(sum(x.^2,2)), Disps,'UniformOutput',false);
    [Tracks.('DistFromWound')] = Dists{:};
    
    RadV = cellfun(@(x) 0.5\diff(x), Dists,'UniformOutput',false);
    [Tracks.('RadialVelocity')] = RadV{:};

    %Velocity - 
    Velocity = arrayfun(@(x) 0.5\([x.tracksCoordAmpCG(9:8:end)',x.tracksCoordAmpCG(10:8:end)',x.tracksCoordAmpCG(11:8:end)']-[x.tracksCoordAmpCG(1:8:end-8)',x.tracksCoordAmpCG(2:8:end-8)',x.tracksCoordAmpCG(3:8:end-8)']),Tracks,'UniformOutput',false);
    SmoothedVelocity = cellfun(@(x) [Smoothing(x(:,1)),Smoothing(x(:,2)),Smoothing(x(:,3))],Velocity,'UniformOutput',false);
    [Tracks.('Velocity')] = SmoothedVelocity{:};
    
     Speed = cellfun(@(x) sqrt(sum(x.^2,2)), SmoothedVelocity,'UniformOutput',false);
     [Tracks.('Speed')] = Speed{:};

    R.setTracks(Tracks,R.PosNames{WellNum})
end
    clearvars Densities Tracks Speed Velocity RadV Dists Disps epiRank CorneaCells;

%%
R.saveResults
%At this point, we should have a results object with PIV data and single
%cell data. Party time.



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


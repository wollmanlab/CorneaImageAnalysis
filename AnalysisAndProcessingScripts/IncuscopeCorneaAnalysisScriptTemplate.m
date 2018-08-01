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
for j=1:numel(R.PosNames);
    PIV = R.getPIVLbl(R.PosNames{j});
    DistanceFromWound = zeros(size(PIV.X,1),numel(PIV.flist));
    RadialVelocity = zeros(size(PIV.X,1),numel(PIV.flist));
    TangentVelocity = zeros(size(PIV.X,1),numel(PIV.flist));
    Speed = zeros(size(PIV.X,1),numel(PIV.flist));
    
    Cent = PIV.WoundLbl.SmoothCentroid(1:numel(PIV.flist)+1);
    %Drift transforms
    Tforms = cell2mat(PIV_MD.getSpecificMetadata('driftTform','Type','EDOF', 'Channel','DeepBlue', 'Position', Wells{j}));
    Tforms = Tforms(:,7:8);
    for i=1:numel(PIV.flist)
        sprintf('frame %d out of %d on position %d out of %d',i,numel(PIV.flist),j,numel(R.PosNames))
        
        x = PIV.X -Cent(i,1);
        y = PIV.Y -Cent(i,2);

        CosAlpha = x./(sqrt(x.^2+y.^2));
        CosAlpha(isnan(CosAlpha))=0;
        SinAlpha = y./(sqrt(x.^2+y.^2));
        SinAlpha(isnan(SinAlpha))=0;
        
        DistanceFromWound(:,i) = sqrt(x.^2+y.^2);
        UMasked = PIV.UMasked(i);%drift correction on u and v
        VMasked = PIV.VMasked(i);
        
        woundU = mean(UMasked);%Tforms(i+1,1)-Tforms(i,1);%Cent(i+1,1)-Cent(i,1);%drift correction on wound motion
        woundV = mean(VMasked);%Tforms(i+1,2)-Tforms(i,2);%Cent(i+1,2)-Cent(i,2);
        
        RadialVelocity(:,i) = (CosAlpha.*(UMasked-woundU)+SinAlpha.*(VMasked-woundV)).*~PIV.indInWound(i);
        TangentVelocity(:,i) = (SinAlpha.*(UMasked-woundU)-CosAlpha.*(VMasked-woundV)).*~PIV.indInWound(i);
        UMasked = UMasked+Tforms(i+1,1)-Tforms(i,1);%drift correction on u and v
        VMasked = VMasked+Tforms(i+1,2)-Tforms(i,2);
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
for WellNum=4%1:numel(Wells)
    pos = Wells{WellNum};
    CorneaCelllbl = cell(numel(frames),1);
    parfor i=R.Frames'
        CorneaCelllbl{i} = CorneaCellsConstructor(fpath, pos,i)
    end
    R.setCorneaCellsLbl(CorneaCelllbl,pos)
end
%% Link adjecent frames
for WellNum=1:numel(Wells)  
    R.Link(R.PosNames(WellNum))
end
%% Close gaps
for WellNum=1:numel(Wells)
    R.closeGaps(R.PosNames(WellNum))
end
%%
R.saveResults
%At this point, we should have a results object with PIV data and single
%cell data. Party time.

%%
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
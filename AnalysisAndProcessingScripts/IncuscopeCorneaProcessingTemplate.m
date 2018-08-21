%% This script loads the raw image stacks, creates SAR, makes drift correction, creates 2d EDoF projections,
%%makes a new MD for these, run a batch PIV analysis for the sequence and
%%do post analysis
BaseStr = regexprep([char(ispc.*'Z:\Images2018\') char(isunix.*'/bigstore/Images2018/')],char(0),'');
Usr = 'Jen';
Project = 'CorneaCCM';
Dataset = 'WoundAgarTitr_2018Mar16_2018May01';
acquisition = 2;
%% Get MD of raw data, open images and resave as SAR
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath);
fpathresults = [fpath filesep 'SAR']
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
Zindexes = max(unique(cell2mat(MD.getSpecificMetadata('Zindex'))));


TSs = MD.getSpecificMetadata('TimestampImage'); %get all image timestamps as unique identifiers
ImgFiles = MD.ImgFiles;
pth = MD.pth;
parfor i=1:size(TSs,1)
    
    foldname = fullfile(pth,'SAR',ImgFiles{i}(1:strfind(ImgFiles{i},'\')-1));
    if ~exist(foldname,'dir')
        mkdir(foldname)
    end
    
    newFilename = fullfile(pth,'SAR',ImgFiles{i}); %generate a new filename!
    newFilename = regexprep(newFilename,'\\',filesep);
    
    %create SAR image
    if  ~exist(newFilename, 'file')
        MD=Metadata(fpath);
        Data = MD.stkread('TimestampImage',TSs{i}, 'flatfieldcorrection', false); %load single image
        wlDecomp = awt2Dlite(Data,4);
        wlDecomp = wlDecomp.*(wlDecomp>0);
        D2 = 0.5*wlDecomp(:,:,:,2)+wlDecomp(:,:,:,3)+0.5*wlDecomp(:,:,:,4);
        imwrite(uint16(D2*2^16),newFilename);
    else
        disp('skipping')
    end
end


%% Drift correction, here, we'll calculate the drift in a few of the planes and just take the mean. 
%If we really want, we could calculate a distinct drift for each frame
%specifically. That's probably unnecessary and computationally expensive.

Channels = {'DeepBlue', 'Red'};
Zindexes = max(unique(cell2mat(MD.getSpecificMetadata('Zindex'))));
ZsToLoad = 1:10:Zindexes-5;
    
for j=1:numel(Wells)
    driftXY.dX = 0;
    driftXY.dY = 0;
    for i=1:numel(frames)-1;
        Data = stkread(MD,'Channel',Channels{1}, 'flatfieldcorrection', false, 'frame', [i i+1], 'Position', Wells{j}, 'Zindex', ZsToLoad);
        imSize = size(Data);
        Data = reshape(Data,[imSize./[1,1,2],2]);
        imSize = size(Data);
        
        dX = [];
        dY = [];
        
        for i=1:imSize(3)
            currRef = Data(:,:,i,1);
            currCorr = Data(:,:,i,2);
            %if you don't value your time, feel free to use imregtform or
            %whatever else you see fit
            imXcorr = convnfft(currRef - mean(currRef(:)),rot90(currCorr,2)-mean(currCorr(:)),'same');
            [maxCorrX, maxCorrY] = find(imXcorr == max(imXcorr(:)));
            dX = [dX maxCorrX-size(currRef,1)/2];
            dY = [dY maxCorrY-size(currRef,2)/2];
        end
        driftXY.dX = [driftXY.dX round(mean(dX))];
        driftXY.dY = [driftXY.dY round(mean(dY))]
    end
    CummulDriftXY.dX = cumsum(driftXY.dX);
    CummulDriftXY.dY = cumsum(driftXY.dY);
    
    
    
    %% Add drift to MD
    Typ = MD.Types;
    Vals = MD.Values;
    
    if ~any(strcmp('driftTform',Typ))
        Typ{end+1}='driftTform'; %Will become a standard in MD.
    end
    Ntypes = size(Typ,2);
    % put the right drift displacements in the right place
    for i=1:numel(frames)
        i
        inds = MD.getIndex({'frame', 'Position'},{i, Wells{j}});
        for j1=1:numel(inds)
            Vals{inds(j1),Ntypes} = [1 0 0 , 0 1 0 , CummulDriftXY.dY(i), CummulDriftXY.dX(i), 1];
        end
    end
    MD.Types = Typ;
    MD.Values = Vals;
end

%% Save new MD in its new home
mdpth = MD.pth;
MD.pth = fullfile(MD.pth,'SAR');
MD.saveMetadataMat;



%% Get MD of SAR data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname filesep 'SAR'];
MD=Metadata(fpath);

fpathresults = [MD.pth filesep 'EDoFandPIVProcessing'];
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
Zindexes = max(unique(cell2mat(MD.getSpecificMetadata('Zindex'))));


%% Run EDOF on each of the hoescht stacks. Save the edof image, the map,
%and the edof of the PI channel using the same map. Also, make a movie of
%the whole timelapse. Do this all in parallel for different positions
%parpool('local', numel(Wells));
EDoFChannel = {'DeepBlue'};
DependentChannels = {'Red'};

parfor j=1:length(Wells)
    
    EDoFChannel = {'DeepBlue'};
    DependentChannels = {'Red'};
    
    MD=Metadata(fpath);%reload metadata in each worker's world so that fileids are right when stkreading
    %determine color ranges
    range = 1:max(cell2mat(MD.getSpecificMetadata('frame')));
    
    
    for i=1:length(range)
        if ~isdir([fpathresults filesep Wells{j}])
            mkdir([fpathresults filesep Wells{j}])
        end
        sprintf('%d of %d', i, length(range))
        Data = stkread(MD,'Channel',EDoFChannel{1}, 'flatfieldcorrection', false, 'frame', i, 'Position', Wells{j});
        Nbins = 2^6; %how many bins to put image in
        if i==1
            [stk, stkmap] = imbedoLMPProp_v3(Data,Nbins); %% For SAR images
        else %move previous map according to registration
            stkmap(stkmap<=5) = Zindexes;
            [stk, stkmap] = imbedo_ByBrightnessMPProp(Data,2^7,stkmap);      
        end
        imwrite(uint16(stk*2^16),sprintf('%sEDOF_%s_%3.3d%s',[fpathresults filesep Wells{j} filesep],EDoFChannel{1},i,'.tiff'));
        imwrite(uint8(255*stkmap/Zindexes),sprintf('%sEDOF_%s_%3.3d%s',[fpathresults filesep Wells{j} filesep],'FocusMap',i,'.tiff'));
        for ChanNum=1:numel(DependentChannels)
            Data = stkread(MD,'Channel',DependentChannels{ChanNum}, 'flatfieldcorrection', false, 'frame', i, 'Position', Wells{j});
            stk = imbedo_ReconstLMP(Data,stkmap);
            imwrite(uint16(stk*2^16),sprintf('%sEDOF_%s_%3.3d%s',[fpathresults filesep Wells{j} filesep],DependentChannels{ChanNum},i,'.tiff'));
        end
        
    end;
end

%% Make Metadata of EDoF images
PIV_MD = Metadata();
PIV_MD.pth = fpathresults;
Channels = [EDoFChannel DependentChannels];
Zindexes = max(unique(cell2mat(MD.getSpecificMetadata('Zindex'))));

for j=1:length(Wells);
    fpthpos = [fpathresults filesep Wells{j}];
    %Channels
    for ChanNum=1:numel(Channels)
        Channel = Channels{ChanNum};
        Fluorophore = unique(MD.getSpecificMetadata('Fluorophore','Channel',Channel));
        Marker = unique(MD.getSpecificMetadata('Marker','Channel',Channel));
        flistDo = getFlistbyPattern(sprintf('%s%s','EDOF_',Channel), fpthpos);
        tf = cell2mat(MD.getSpecificMetadata('TimestampFrame','Position',Wells{j},'Channel',Channel));
        tf = mean(reshape(tf,[],numel(flistDo)));
        driftTforms = MD.getSpecificMetadata('driftTform','Position',Wells{j},'Channel','DeepBlue', 'Zindex', 1);
        PixelSize = MD.getSpecificMetadata('PixelSize','Position',Wells{j},'Channel','DeepBlue', 'Zindex', 1);
        
        for i=1:length(flistDo)
            PIV_MD.addNewImage([Wells{j} filesep flistDo{i}],'TimestampFrame',tf(i),'Type','EDOF','Channel',Channel,'Fluorophore',Fluorophore,'Marker',Marker,'Position', Wells{j},'frame',i,'PixelSize',PixelSize, 'driftTform',driftTforms{i},'Zplanes',Zindexes);
        end;
    end
    
    %Maps
    Channel = 'All';
    Fluorophore = 'All';
    Marker = 'All';
    flistDo = getFlistbyPattern('FocusMap', fpthpos);
    tf = cell2mat(MD.getSpecificMetadata('TimestampFrame','Position',Wells{j},'Channel',Channels{1}));
    tf = mean(reshape(tf,[],numel(flistDo)));
    driftTforms = MD.getSpecificMetadata('driftTform','Position',Wells{j},'Channel','DeepBlue', 'Zindex', 1);
    PixelSize = MD.getSpecificMetadata('PixelSize','Position',Wells{j},'Channel','DeepBlue', 'Zindex', 1);
    
    for i=1:length(flistDo)
        PIV_MD.addNewImage([Wells{j} filesep flistDo{i}],'TimestampFrame',tf(i),'Type','FocusMap','Channel',Channel,'Fluorophore',Fluorophore,'Marker',Marker,'Position', Wells{j},'frame',i,'PixelSize',PixelSize, 'driftTform',driftTforms{i},'Zplanes',Zindexes);
    end;   
end;

PIV_MD.saveMetadataMat(fpathresults);

%% save registered pics for PIV
EDoFChannel = {'DeepBlue'};

for j=1:numel(Wells)
    Data = stkread(PIV_MD,'Channel',EDoFChannel{1},'Type', 'EDOF', 'flatfieldcorrection', false, 'Position', Wells{j},'register',true);
    for i=1:size(Data,3)
        imwrite(uint16(Data(:,:,i)*2^16),sprintf('%sReg_EDOF_%s_%3.3d%s',[fpathresults filesep Wells{j} filesep],EDoFChannel{1},i,'.tiff'));
    end
end

%% Get MD of SAR data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname filesep 'SAR'];
MD=Metadata(fpath);

fpathresults = [MD.pth filesep 'EDoFandPIVProcessing'];
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
Zindexes = max(unique(cell2mat(MD.getSpecificMetadata('Zindex'))));
PixelSize = MD.getSpecificMetadata('PixelSize','Channel','DeepBlue', 'Zindex', 1);
PixelSize = PixelSize{1};


%% Load MD of EDoF images
fpathresults = [fpath filesep 'EDoFandPIVProcessing']
PIV_MD = Metadata(fpathresults);
Wells = unique(PIV_MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(PIV_MD.getSpecificMetadata('frame')));




%% New way to make movies
j=1
Data = stkread(PIV_MD,'sortby','Channel','Type', 'EDOF', 'flatfieldcorrection', false, 'Position', Wells{j},'resize', 0.5,'register',true);

stkshow(Data);
Channels = unique(PIV_MD,'Channel');
PixelSize = PIV_MD.getSpecificMetadataByIndex('PixelSize', 1);
PixelSize = PixelSize{1}{1};

MIJ.selectWindow('Data');
MIJ.run('Stack Splitter', 'number=2');%Split stack to 2 colors
MIJ.selectWindow('Data');
MIJ.selectWindow('stk_0002_Data');
MIJ.run('Enhance Contrast...', 'saturated=0.2 normalize process_all use');%stretch to range and fix LUT
MIJ.selectWindow('stk_0001_Data');
MIJ.run('Enhance Contrast...', 'saturated=0.2 normalize process_all use');
MIJ.run('Concatenate...', 'image1=stk_0001_Data image2=stk_0002_Data');%merge colors back
MIJ.run('8-bit');%convert to 8 bit
MIJ.run('Stack Contrast Adjustment', 'is');%equalize contrast between frames
MIJ.run('Stack to Hyperstack...', ['order=xytcz channels=' num2str(numel(Channels)) ' slices=1 frames=' num2str(size(Data,3)/numel(Channels)) ' display=Composite']); %Make into xytc hyperstack
MIJ.run('Properties...', ['channels=' num2str(numel(Channels)) ' slices=1 frames=' num2str(size(Data,3)/numel(Channels)) ' unit=um pixel_width=' num2str(PixelSize) ' pixel_height=' num2str(PixelSize) ' voxel_depth=1.0000']);%set pixel size
MIJ.run('Scale Bar...', 'width=100 height=8 font=28 color=White background=None location=[Lower Right] bold overlay');%add scale bar
MIJ.run('Time Stamper', ['starting=0 interval=0.5 x=0 y=50 font=50 decimal=1 anti-aliased or=h overlay']);%add time stamp

f = figure;%wait for user to adjust colors, intensity, etc
set(f,'Name','Please adjust brightness and contrast','Position',[360 638 240 60],'NumberTitle', 'off', 'Toolbar','none','Menubar','none')
h = uicontrol('Position',[20 20 200 40],'String','OK',...
              'Callback','uiresume(gcbf)');
uiwait(f); 
close(f);
%save
MIJ.run('AVI... ', ['compression=JPEG frame=14 save=' sprintf('%sMovie_EDOF_%s%s',[fpathresults filesep Wells{j} filesep],Wells{j},'.avi')]);

 
%% Load MD of EDoF images
fpathresults = [fpath filesep 'EDoFandPIVProcessing']
PIV_MD = Metadata(fpathresults);
Wells = unique(PIV_MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(PIV_MD.getSpecificMetadata('frame')));
%% Open jpiv and make sure all parameters are right
    fpthpos = [fpathresults filesep Wells{1}];

    system(sprintf('cd %s; java -jar -Xmx16G /home/alon/bin/JPIV/jpivlib/jpivc.jar %s &',fpthpos))
    %General: Consecutive 1+2 2+3 ...
    %Window:
    %Multi pass: 3
    %Width\Height: 128, 64, 64
    %Search width\height: 64, 48, 32
    %Spacing: 128, 64, 64
    %Normalize median test + replace
%% Run batch PIV analysis by sequence ([1-2], [2-3]...)
for j=1:length(Wells)
    fpthpos = [fpathresults filesep Wells{j}];
    flistDo = getFlistbyPattern('Reg_EDOF_DeepBlue', fpthpos);
    system(sprintf('cd %s; java -jar -Xmx16G /home/alon/bin/JPIV/jpivlib/jpivc.jar %s &',fpthpos,strjoin(flistDo,' ')))
end;
%% Do post-analysis of PIV - removing spurious vectors, etc.
%results will be saved in the same folder as PIV_Results_001,2,3...
for j=1:length(Wells)
    fpthpos = [fpathresults filesep Wells{j} filesep];
    postAnalysis(fpthpos)
end;

%end of Processing Script





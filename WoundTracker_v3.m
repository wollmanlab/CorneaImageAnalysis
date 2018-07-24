function Woundlbl = WoundTracker_v3(pth, Well)
close all;

MD=Metadata(pth); 
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));


NucFrame=stkread(MD,'Channel','DeepBlue','Type', 'EDOF','Position',Well,'specific',1,'flatfieldcorrection',false);
%N = 100;%point for wound tracking;
Woundlbl = WoundLbl;
Woundlbl.PosName = Well;
Woundlbl.pth = pth;
Woundlbl.ImageDims = [size(NucFrame,1),size(NucFrame,2)];
Woundlbl.PolyXY =  cell(1,length(frames));
Woundlbl.IsThereAWound = zeros(1,length(frames));
%Woundlbl.Tvec= MD.getSpecificMetadata('TimestampFrame','Position',Well,'Channel','DeepBlue');
w=1;
for i=1:length(frames)-1
        NucFrame=stkread(MD,'Channel','DeepBlue','Position',Well,'specific',i,'flatfieldcorrection',false);
        sprintf('%d out of %d',i,length(frames))
        if w==1;
        imagesc(NucFrame); shg;
        set(gcf,'Position',[10 100 600 480])
            IsThereAWound = questdlg('Is this cornea wounded?', ...
            'Wound?','Yes','No','No');
        end
        
        switch IsThereAWound
            case 'Yes'               
                XY = impoly;
                position = XY.getPosition;
            case 'No'
                w=0;
                continue;     
           close all;
        end
        
   
Woundlbl.PolyXY{i} =  position;
Woundlbl.IsThereAWound(i) = w;
end;

end
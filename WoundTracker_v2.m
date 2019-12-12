function Woundlbl = WoundTracker_v2(pth, Well,varargin)
close all;
channel = ParseInputs('channel','DeepBlue',varargin);

MD=Metadata(pth); 
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));

NucFrame=stkread(MD,'Channel',channel,'Position',Well,'specific',1,'flatfieldcorrection',false);
N = 100;%point for wound tracking;
Woundlbl = WoundLbl;
Woundlbl.PosName = Well;
Woundlbl.pth = pth;
Woundlbl.ImageDims = [size(NucFrame,1),size(NucFrame,2)];
Woundlbl.PolyXY =  zeros(N,2,length(frames));
Woundlbl.IsThereAWound = zeros(1,length(frames));
%Woundlbl.Tvec= MD.getSpecificMetadata('TimestampFrame','Position',Well,'Channel','DeepBlue');

for i=1:length(frames)
sprintf('%d out of %d',i,length(frames))
    if i==1 %initial wound
        imagesc(NucFrame); shg;
        set(gcf,'Position',[10 100 600 480])

        IsThereAWound = questdlg('Is this cornea wounded?', ...
            'Wound?','Yes','No','No');
        switch IsThereAWound
            case 'Yes'
                Cent=ceil(ginput(1));
                w = 1;
            case 'No'
                continue;
        end
        
    end
    close all;
    switch IsThereAWound
        case 'Yes'
            try
                NucFrame=stkread(MD,'Channel',channel,'Position',Well,'specific',i,'flatfieldcorrection',false); 
                if i==1
                    XY = FindWoundPoly(NucFrame,Cent,N,'Present');
                    Circ = sum(sqrt(diff(XY(:,1)).^2+diff(XY(:,2)).^2));
                    Circ0 = Circ;
                else
                    XY = FindWoundPoly(NucFrame,Cent,N);
                    Circ = sum(sqrt(diff(XY(:,1)).^2+diff(XY(:,2)).^2));
                end
            catch 
                Circ = 0;%if there is an error, ask what to do
            end
       eflag = Circ<0.9*Circ0 || Circ>1.2*Circ0;  
       while eflag;%find bad segmentation and ask for human help
           imagesc(NucFrame); shg;
           set(gcf,'Position',[10 100 600 480])
           IsThereAWound = questdlg('We seem to be experiencing some technical difficulties. Is this cornea still wounded?', ...
            'Wound?','Yes','No','No');
            switch IsThereAWound
                case 'Yes'
                    Cent=ceil(ginput(1));
                    XY = FindWoundPoly(NucFrame,Cent,N,'Present');
                    reWound = questdlg('Did it work?','Worked?','Yes','No','No');
                    if strcmp(reWound,'Yes');
                            eflag = 0;
                    end
                case 'No'
                    eflag = 0;
                    continue
            end
        end
                
        case 'No'
            continue
    end
    
Circ0 = Circ;    
    
Woundlbl.PolyXY(:,:,i) =  XY;
Woundlbl.IsThereAWound(i) = w;
end;

end
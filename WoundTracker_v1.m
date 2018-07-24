function Wound = WoundTracker_v1(pth, Well)
close all;



MD=Metadata(pth); 
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));

NucFrame=stkread(MD,'Channel','DeepBlue','Position',Well,'specific',frames(1),'flatfieldcorrection',false);

Wound = struct;
Wound.Masks =  ones(size(NucFrame,1),size(NucFrame,2),length(frames));
Wound.Centroid = repmat(ceil(size(NucFrame)./2),length(frames),1);
Wound.Area = zeros(1,length(frames));
Wound.IsThereAWound = zeros(1,length(frames));


for i=1:length(frames)
sprintf('%d out of %d',i,length(frames))
    if i==1 %initial wound
        imagesc(NucFrame); shg;
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
                NucFrame=stkread(MD,'Channel','DeepBlue','Position',Well,'specific',frames(i),'flatfieldcorrection',false); 
                [bw,Cent,Area] = FindWound(NucFrame,Cent);
                if i==1
                    Area0 = Area;
                end
            catch 
                Area = 0;%if there is an error, ask what to do
            end
       eflag = Area<0.75*Area0 || Area>1.5*Area0;  
       while eflag;%find bad segmentation and ask for human help
           imagesc(NucFrame); shg;
           IsThereAWound = questdlg('We seem to be experiencing some technical difficulties. Is this cornea still wounded?', ...
            'Wound?','Yes','No','No');
            switch IsThereAWound
                case 'Yes'
                    Cent=ceil(ginput(1));
                    [bw,Cent,Area] = FindWound(NucFrame,Cent,'Present');
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
    
Area0 = Area;    
    
Wound.Masks(:,:,i) =  ~bw;
Wound.Centroid(i,:) = Cent;
Wound.Area(i) = Area;
Wound.IsThereAWound(i) = w;
end;

end
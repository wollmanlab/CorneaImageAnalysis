classdef CorneaCellsLbl < handle %class of single cell processing of a whole cornea at a single timepoint
    properties
        PosName
        Frame
        pth
        ImageDims
        zAspect
        
        Centroids
        Intensities
        num
        
        TopoZ
        %GridZ 
        CC

        Link12Mat
        Link21Mat
    end
    
    properties (Transient = true)
        verbose = true;
    end
    
    properties (Dependent = true)
        filenames
    end
    
    
    
    methods
        function scatter(W,varargin)
            tzeva = viridis(length(W.Centroids));
            if nargin==1
                scatter3(W.Centroids(:,1),W.Centroids(:,2),-W.Centroids(:,3),[],tzeva);
            else
                range = varargin{1};
                if any(range>W.num)
                    range = range(1):W.num;
                    disp('range out of bounds, drawing all points up to the total # of cells.')
                end;
                scatter3(W.Centroids(range,1),W.Centroids(range,2),-W.Centroids(range,3),[],tzeva(range,:));
            end
            %scatter(Centroids(:,1),Centroids(:,2),[],parula(length(Centroids)));
            set(gca,'xlim',[-200 3000],'ylim',[-200 3000],'CameraPositionMode','manual','CameraPosition',[-5685       -6469        3001])
            %hold on
            shg
            drawnow
            pause(0.01)
        end
        
        function stkshow(W)
            MD=Metadata(W.pth);
            
            Data =  stkread(MD,'Channel','DeepBlue', 'flatfieldcorrection', false, 'frame', W.Frame, 'Position', W.PosName,'register',false);
            RChannel=Data;
            GChannel=Data;
            BChannel=Data;
            [h,x] = hist(log(datasample(Data(:),min(1000000,numel(Data(:))))),1000);
            maxC = exp(x(find(cumsum(h)./sum(h)>0.995,1,'first')));
            cmap = parula(W.CC.NumObjects)*maxC*1.5;%scale colormap so that data and centroids are all visible
            
            for i=1:W.CC.NumObjects
                indexToChange = W.CC.PixelIdxList{i};
                RChannel(indexToChange)=cmap(i,1); %Replace actual pixel values with relevant RGB value
                GChannel(indexToChange)=cmap(i,2);
                BChannel(indexToChange)=cmap(i,3);
                i
            end
            RGB = cat(3,RChannel,GChannel,BChannel);%combine into a single RGB image and show. stkshow is sloooooooooooowing me down.
            stkshow(RGB);
            MIJ.selectWindow('RGB');
            MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=3 slices=' num2str(size(Data,3)) ' frames=1 display=Composite']);
            %stkshow(RChannel)
            %stkshow(GChannel)
            %stkshow(BChannel)
            %MIJ.run('Merge Channels...', 'c1=RChannel c2=GChannel c3=BChannel create');            
        end
        
        
        function scattershow(W)
            MD=Metadata(W.pth);
            
            Data =  stkread(MD,'Channel','DeepBlue', 'flatfieldcorrection', false, 'frame', W.Frame, 'Position', W.PosName,'register',false);
            RChannel=zeros(size(Data));
            GChannel=zeros(size(Data));
            BChannel=zeros(size(Data));
            [h,x] = hist(log(datasample(Data(:),min(1000000,numel(Data(:))))),1000);
            maxC = exp(x(find(cumsum(h)./sum(h)>0.995,1,'first')));
            cmap = parula(W.CC.NumObjects)*maxC*1.5;%scale colormap so that data and centroids are all visible
            
            for i=1:W.CC.NumObjects
                indexToChange = W.CC.PixelIdxList{i};
                RChannel(indexToChange)=cmap(i,1); %Replace actual pixel values with relevant RGB value
                GChannel(indexToChange)=cmap(i,2);
                BChannel(indexToChange)=cmap(i,3);
                i
            end
            RGB = cat(3,RChannel,GChannel,BChannel);%combine into a single RGB image and show. stkshow is sloooooooooooowing me down.
            stkshow(RGB);
            MIJ.selectWindow('RGB');
            MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=3 slices=' num2str(size(Data,3)) ' frames=1 display=Composite']);
            %stkshow(RChannel)
            %stkshow(GChannel)
            %stkshow(BChannel)
            %MIJ.run('Merge Channels...', 'c1=RChannel c2=GChannel c3=BChannel create');            
        end
        
    end
end

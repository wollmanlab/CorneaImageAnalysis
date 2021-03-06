function VSmooth = Smoothing(V,varargin)
    %V = V-mean(V);
    shp = size(V);
    if ~logical(length(shp) == 2 && min(shp) == 1);
        VSmooth = zeros(shp);
        for i=1:shp(2)
            VSmooth(:,i) = Smoothing(V(:,i));
        end
        %error('Input must be a vector')
    end
    V = V(:);
    
    
    nanInds = find(isnan(V));
    if any(nanInds)%fill gaps with linear interpolation and dashed lines
        %warning('Filling NaNs with linear intepolation')

        skipInds = cumsum([1; diff(find(isnan(V)))>1]);
        
        for ind=1:numel(unique(skipInds))
            jx = skipInds==ind;
            indsToFill = nanInds(jx);
            if indsToFill(1)==1 %deal with edge cases by nn interp
                V(indsToFill(1):indsToFill(end)+1) = V(indsToFill(end)+1);
            elseif indsToFill(end)==length(V)
                V(indsToFill(1)-1:indsToFill(end)) = V(indsToFill(1)-1);
            else
            m = (V(indsToFill(end)+1)-V(indsToFill(1)-1))/(1+sum(jx));
            V(indsToFill(1)-1:indsToFill(end)+1) = V(indsToFill(1)-1)+m*(0:sum(jx)+1);
            end
        end

    end

    

    
    arg.neigh = 5;
    
    arg = parseVarargin(varargin,arg);
    
    
    frame = 1:length(V);
    
    Neigh = arg.neigh;
    
    flt = GaussianFit([1, 0, Neigh/2], -Neigh:Neigh);
    Cent = zeros(numel(frame),2);
    PadBefore = -(frame(1)-Neigh-1)*((frame(1)-Neigh)<1);
    PadAfter =(frame(end)+Neigh-length(V))*((frame(end)+Neigh)>length(V));
    CentroidEnv = [repmat(V(1),PadBefore,1);...
        V(max(frame(1)-Neigh,1):min(frame(end)+Neigh,length(V)));...
        repmat(V(length(V)),PadAfter,1)];
    VSmooth = [];
    for i=1:numel(frame)
        VSmooth = [VSmooth, sum(flt'.*CentroidEnv(i:i+Neigh*2,:))];
    end
    VSmooth = reshape(VSmooth,shp);
    %plot(1:240,VSmooth, 1:240, V); shg
end

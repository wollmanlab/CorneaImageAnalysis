function VSmooth = Smoothing(V,varargin)
    %V = V-mean(V);
    
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
    
    %plot(1:240,VSmooth, 1:240, V); shg
end

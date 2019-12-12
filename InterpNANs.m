function V = InterpNANs(V,varargin)
    %V = V-mean(V);
    shp = size(V);
    if ~logical(length(shp) == 2 && min(shp) == 1);
        error('Input must be a vector')
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

end

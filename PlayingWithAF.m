 flt = fspecial('average');
 U = filter2(flt,fmap1);
 [h,x] = histcounts(reshape((sqrt((fmap1-U).^2)),[],1),100);
 J = cumsum(h)/sum(h)>.95;
 ind2switch = find(sqrt((fmap1-U).^2) >= min(x(J)));
 %%
 for repInd = 1:size(ind2switch)
     [I,J] = ind2sub([size(fmap1,1) size(fmap1,2)], ind2switch(repInd));
    if I>1 & I<size(fmap1,1) & J>1 & J<size(fmap1,2)
     fmap1(I,J) = (fmap1(I+1,J)+fmap1(I-1,J)+fmap1(I,J+1)+fmap1(I,J-1)+...
                      fmap1(I+1,J+1)+fmap1(I-1,J+1)+fmap1(I+1,J-1)+fmap1(I-1,J-1))/8;
                  
                  %sides
    elseif I==1 & J>1 & J<size(fmap1,1)
     fmap1(I,J) = (fmap1(I+1,J)+fmap1(I,J+1)+fmap1(I,J-1)+...
                      fmap1(I+1,J+1)+fmap1(I+1,J-1))/5;
    
    elseif I==size(fmap1,1) & J>1 & J<size(fmap1,2)
     fmap1(I,J) = (fmap1(I-1,J)+fmap1(I,J+1)+fmap1(I,J-1)+...
                      fmap1(I-1,J+1)+fmap1(I-1,J-1))/5;        
    
    elseif I>1 & I<size(fmap1,1) & J==1
     fmap1(I,J) = (fmap1(I+1,J)+fmap1(I-1,J)+fmap1(I,J+1)+...
                      fmap1(I+1,J+1)+fmap1(I-1,J+1))/5;
    
    elseif I>1 & I<size(fmap1,1) &  J==size(fmap1,2)
     fmap1(I,J) = (fmap1(I+1,J)+fmap1(I-1,J)+fmap1(I,J-1)+...
                      fmap1(I+1,J-1)+fmap1(I-1,J-1))/5;
                  %corners
    elseif I==1 & J==1
     fmap1(I,J) = (fmap1(I+1,J)+fmap1(I,J+1)+...
                      fmap1(I+1,J+1))/3;
    
    elseif I==size(fmap1,1) & J==1
     fmap1(I,J) = (fmap1(I-1,J)+fmap1(I,J+1)+...
                      fmap1(I-1,J+1))/3;        
    
    elseif I==1 &  J==size(fmap1,2)
     fmap1(I,J) = (fmap1(I+1,J)+fmap1(I,J-1)+...
                      fmap1(I+1,J-1))/3;
    
    elseif I==size(fmap1,1) &  J==size(fmap1,2)
     fmap1(I,J) = (fmap1(I-1,J)+fmap1(I,J-1)+...
                      +fmap1(I-1,J-1))/3;          
    end;  
 end

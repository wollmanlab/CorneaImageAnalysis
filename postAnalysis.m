function postAnalysis(fpath)
    cd(fpath)
    Temp = 'deleteme';
    endTemp = '.jvc';
    endTempNo = '.jvc~';
    flist = dir(fpath);
    flist = {flist.name};

    ix1 = regexp(flist, endTemp);
    ix2 = regexp(flist, Temp);
    ix3 = regexp(flist, endTempNo);
    ix=~cellfun('isempty',ix1) & ~cellfun('isempty',ix2) & cellfun('isempty',ix3) ;
    flistDo = flist(ix);
    flistDo'
    
    flt = fspecial('disk',1);
    %flt = fspecial('average');
    %flt = [0 1 0; 1 1 1; 0 1 0]./5;
%%
for ii=1:length(flistDo) 
    sprintf('%d of %d', ii, length(flistDo))
    
    data = load(flistDo{ii});
    x = data(:,1);
    y = data(:,2);
    u = data(:,3);
    v = data(:,4);
    s = data(:,5);
    delete(flistDo{ii});
    %% Get number of vectors
    % Get number of x and y vectors - searches through y column for a change in number. 
    % This gives number of x co-ordinates, nox. When the number changes, the
    % difference, k,  is found. Noy is found by subtracting the co-ordinate of the
    % first vector from the last, then dividing by k and adding 1

    noy = sum(diff(y)~=0)+1;
    nox=size(y,1)/noy;
    %% Convert Vectors to Matrices
    U=reshape(u,nox,noy);
    V=reshape(v,nox,noy);
    X=reshape(x,nox,noy);
    Y=reshape(y,nox,noy);
    S=reshape(s,nox,noy);
    %% filter noisy vectors
    %U = imfilter(U,flt,'replicate');
    %V = imfilter(V,flt,'replicate');
    
%OK, interpolation introduces systematic errors. in stead, we will just remove
%the 2% of vectors that are least like their neighbors.
 flt = [1 1 1; 1 0 1; 1 1 1]./8;

 Uneigh = imfilter(U,flt,'replicate');
 Vneigh = imfilter(V,flt,'replicate');
 [h,x] = histcounts(reshape(sqrt((U-Uneigh).^2+(V-Vneigh).^2),[],1),100);
 J = cumsum(h)/sum(h)>.98;
 ind2switch = find(sqrt((U-Uneigh).^2+(V-Vneigh).^2) >= min(x(J)));


    %% return to vector form and save
    x = reshape(X,nox*noy,1);
    y = reshape(Y,nox*noy,1);
    u = reshape(U,nox*noy,1);
    v = reshape(V,nox*noy,1);
    s = reshape(S,nox*noy,1);

    u(ind2switch) = 0;
    v(ind2switch) = 0;
    % write output file
    dataOut = [x, y, u, v, s];
    dlmwrite(sprintf('%sPIV_Results_%3.3d.jvc',fpath, ii), dataOut, 'delimiter', ' ','precision', '%+.4E');
end
end
 function XY = FindWoundPoly(NucFrame,xystart,N,varargin)

j1=find(strcmp(varargin, 'Present'));
if length(j1)>1,
    error('Too many present inputs!')
elseif length(j1)==1,
    Present = true;
else
    Present = false;
end;

 %%
 
%N=100; 
Iter=100; 
xystart=ceil(xystart); 

R=zeros(size(NucFrame)); 
R(xystart(2),xystart(1))=1; 
R=bwdist(R); 

[Xout,Yout]=meshgrid((1:size(NucFrame,2))-xystart(1),(1:size(NucFrame,1))-xystart(2));
Xout=Xout./R; 
Yout=Yout./R;
Xout(xystart(2),xystart(1))=0; 
Yout(xystart(2),xystart(1))=0; 

NucMap = 1-mat2gray(imfilter(NucFrame,fspecial('gauss',150,25),'replicate')); 
r=15; 
a=linspace(-pi,pi,N+1)';
a(end)=[]; 


NucFlt = imfilter(NucMap,fspecial('gauss',150,20),'replicate');
[NucFlt_dX,NucFlt_dY]=gradient(NucFlt+8*exp(-R/5)); %Added some regulation


%%
XY=repmat(xystart,N,1)+r.*[sin(a) cos(a)]; 

wF1 = 1.2; 
wF2 = -0.8; 
wF3 = -1200; 

Circ = 0;
%%
for i=1:3000
    %%
    
    ix=sub2ind(size(NucMap),ceil(XY(:,2)),ceil(XY(:,1))); 
%     Vnuc=NucMap(ix); 
    
    
    for ii=1:N, 
        if ii==1,
            abv=N; 
        else
            abv=ii-1; 
        end, 
        if ii==N,
            blw=1;
        else
            blw=ii+1;
        end
        XYn(ii,:)=mean(XY([abv blw],:));
    end
    F1=[Xout(ix) Yout(ix)];
    F2=(XY-XYn);    
    F3 = [NucFlt_dX(ix) NucFlt_dY(ix)]; 
    
    Ftot = wF1*F1+wF2*F2+wF3*F3;
    
    Circ0 = Circ; %Using circumference for stopping criteria
    XY=XY+Ftot; 
    
    Circ = sum(sqrt(diff(XY(:,1)).^2+diff(XY(:,2)).^2));

    if abs((Circ-Circ0)./Circ)<5E-5 || abs(Circ)<150;
        if Present
            figure(1)
            hold on
            plot(XY(:,1),XY(:,2),'.-') 
        end
        break
    end
    %% from here on is just plotting
    if Present
        if mod(i,50)
            continue
        end
    dF1=wF1*sqrt(sum(F1.^2,2)); 
    dF2=wF2*sqrt(sum(F2.^2,2)); 
    dF3=wF3*sqrt(sum(F3.^2,2)); 
    dF=sqrt(sum(Ftot.^2,2)); 
    figure(1)
    set(gcf,'Position',[10 100 600 480])
    clf
    imagesc(NucMap); 
    % axis xy
    hold on
    plot(XY(:,1),XY(:,2),'.-')
    
    abs((Circ-Circ0)./Circ)
    
    figure(2)
    set(gcf,'Position',[620 200 600 240])
    clf
    subplot(1,2,1)
    hold on
    plot(XY(:,1),XY(:,2),'.',XYn(:,1),XYn(:,2),'x')
    quiver(XY(:,1),XY(:,2),F1(:,1)*20,F1(:,2)*20,'b')
    quiver(XY(:,1),XY(:,2),F2(:,1)*5,F2(:,2)*5,'r')
    axis([0 2448 0 2048])
    axis equal
    title(num2str(i))
    
    subplot(2,2,2)
    bar([dF1 dF2 dF3],'stacked')
    ylim([0 1])
    xlim([0, N])
    subplot(2,2,4)
    bar(dF)
    ylim([0 1])
    xlim([0, N])
    drawnow
    end
end

%bw = poly2mask(double(XY(:,1)),double(XY(:,2)),size(NucFrame,1),size(NucFrame,2));
%Wound = regionprops(bw,'Centroid','Area');
%Centroid = Wound.Centroid;
%Area = Wound.Area;
 end
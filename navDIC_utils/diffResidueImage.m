global hd
clearvars -except hd

IMG = cat(4,hd.Images{1}{:}) ;
Seed = hd.Seeds(4) ;
Nodes = Seed.MovingPoints(:,:,1) ;

[nI,nJ,nC,nFr] = size(IMG) ;
[JJ,II] = meshgrid(1:nJ,1:nI) ;% Get the Seed

% Get Infos
    Elems = Seed.Triangles ;
    nNodesByElements = size(Elems,2) ;
    nElems = size(Elems,1) ;
    nNodes = size(Nodes,1) ;
    
% Processes
    globalDIC_02_1_ProcessGeometries ;
    globalDIC_02_2_ShapeFunctions ;
    
%%
convFilt = @(I)conv2(I,ones(1),'same') ;
UU = @(ii,comp)MAPPING*(Seed.MovingPoints(:,comp,ii)-Seed.MovingPoints(:,comp,1)) ;
currentImg  = @(ii)reshape(DOMAIN'*interp2(JJ,II,double(IMG(:,:,:,ii)),JJ(indDOMAIN)+UU(ii,1),II(indDOMAIN)+UU(ii,2),'linear',0),[nI nJ]) ;
%diffImg = @(ii)currentImg(ii) ;
diffImg = @(ii)abs(currentImg(ii)-currentImg(ii-1)) ;
%diffImg = @(ii)abs(currentImg(ii)-currentImg(1)) ;

clf
    im = imagesc(IMG(:,:,:,1)) ;
        axis ij
        axis equal
        axis tight
        set(gca,'xtick',[],'ytick',[])
        colorbar
    slider = uicontrol(gcf ...
                ,'style','slider' ...
                ,'Units','normalized' ...
                ,'position',[0 0 1 0.025] ...
                ,'min',2 ...
                ,'max',size(IMG,4) ...
                ,'SliderStep',[1/(size(IMG,4)-1) 1/10] ...
                ...,'callback',@(src,evt)set(im,'CData',diffImg(round(src.Value))) ...
                ) ;
    addlistener(slider,'Value','PostSet',@(src,evt)set(im,'CData',diffImg(round(slider.Value)))) ;
    slider.Value = 2 ;
    
    
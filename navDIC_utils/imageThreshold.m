%% THRESHOLD IMAGES
    global hd
    uiSize = 20 ; 
    margin = 0.001 ; 
    alpha = 0.5 ;
    maxMedFiltRadius = 20 ;
    
    prev = navDICCameraPreview(hd) ;
    hd.Previews(end+1) = prev ;
    
    % Increase figure size
    prev.fig.Position = prev.fig.Position + 2*[-1 0 1 0]*uiSize ;
    
    % Add uicontrols on top of figure
    figSize = prev.fig.Position(3:4) ;
    sliderWidth = uiSize/figSize(1) ; 
    btnHeight = uiSize/figSize(2) ; 
    maxT = max(getrangefromclass(hd.Images{prev.CameraID}{end})) ;
    btnInv = uicontrol(prev.fig,'style','toggle'...
                        ,'Units','normalized' ...
                        ,'Position',[0 1 0 0] + [0 0 1 0]*sliderWidth + margin*[1 1 -2 -2] + [0 -1 0 1]*btnHeight ...
                        ,'Value',0 ...
                        ,'String','~' ...
                        ,'Tooltip','Invert mask' ...
                        ) ;
    tSlider = uicontrol(prev.fig,'style','slider'...
                        ,'Units','normalized' ...
                        ,'Position',[0 0 0 1] + [0 0 1 0]*sliderWidth + margin*[1 1 -2 -2] - [0 0 0 1]*btnHeight ...
                        ,'Min',0,'Max',maxT,'Value',0 ...
                        ,'Tooltip','Threshold Value' ...
                        ) ;
    rSlider = uicontrol(prev.fig,'style','slider'...
                        ,'Units','normalized' ...
                        ,'Position',[0 0 0 .5] + [1 0 1 0]*sliderWidth + margin*[1 1 -2 -2] ...
                        ,'Min',0,'Max',maxMedFiltRadius,'Value',0 ...
                        ,'Tooltip','Stat. Filter Radius' ...
                        ) ;
    oSlider = uicontrol(prev.fig,'style','slider'...
                        ,'Units','normalized' ...
                        ,'Position',[0 .5 0 .5] + [1 0 1 0]*sliderWidth + margin*[1 1 -2 -2] ...
                        ,'Min',0,'Max',1,'Value',.5 ...
                        ,'Tooltip','Stat. Filter Order' ...
                        ) ;
                    
    % Shift the preview axes
    prev.Img.Parent.OuterPosition = [0 0 1 1] + 2*[1 0 -1 0]*sliderWidth ;
    
    % Init mask
    maskIm = copyobj(prev.Img,prev.Img.Parent) ;
    maskIm.CData = maskIm.CData*0 + reshape([1 0 0],[1 1 3]) ;
    
    % Callback definition
    roundMask = @(R)((-R:R).^2 + (-R:R)'.^2)<=R^2 ;
    filtImg = @(img,R,ord)ordfilt2(img,max(round(pi*R^2*ord),1),roundMask(R),'symmetric') ;
    getMask = @(img)filtImg(abs((img>=tSlider.Value)-btnInv.Value),ceil(rSlider.Value),oSlider.Value) ;
    setMask = @()set(maskIm,'AlphaData',getMask(prev.Img.CData(:,:,1))*alpha) ;
    setMask() ;
    
    % Listeners
    btnInv.Callback = @(src,evt)setMask() ;
    addlistener(tSlider,'Value','PostSet',@(src,evt)setMask()) ;
    addlistener(rSlider,'Value','PostSet',@(src,evt)setMask()) ;
    addlistener(oSlider,'Value','PostSet',@(src,evt)setMask()) ;
    addlistener(prev.Img,'CData','PostSet',@(src,evt)setMask()) ;
    prev.updatePreview(hd) ;
    
%% APPLY TO ALL IMAGES
    cam = prev.CameraID ;
    IMG = hd.Images{cam} ;
    wtbr = waitbar(0,'Thresholding...') ;
    for fr = 1:numel(IMG)
        IMG{fr} = getMask(IMG{fr}) ;
        wtbr = waitbar(fr/numel(IMG),wtbr) ;
    end
    delete(wtbr) ;
    
    
%% PUSH TO NAVDIC
    newCamIdx = numel(hd.Cameras)+1 ; 2 ;
    newCam = hd.Cameras(cam) ;
    newCam.Name = ['Threshold | ' hd.Cameras(cam).Name] ;
    hd.Images{newCamIdx} = IMG ;
    hd.Cameras(newCamIdx) = newCam ;
    
    
    
    
    
    
    
    
%% TRANSFORM IMAGES
    global hd
    cam = 2 ;
    refFrame = hd.nFrames ;
    elmt =  ... pkg.geometry.mesh.elements.base.Quadrangle ...
            pkg.geometry.mesh.elements.Quad8 ...
            ;   
    pos = [1 1 ; nJ 1 ; nJ nI ; 1 nI] ;
    gridDiv = 10 ;
        
    refImg = hd.Images{cam}{refFrame} ;
    [nI,nJ,nC] = size(refImg) ;
    [JJ,II] = meshgrid(0:nJ,0:nI) ;
    [GX,GY] = meshgrid(linspace(0,1,gridDiv),linspace(0,1,gridDiv)) ;
        
    clf ;
    ax1 = axes('OuterPosition',[0 0 .5 1]) ;
        imagesc(refImg)
        p1 = (1:gridDiv-1)'+(0:gridDiv-2)*gridDiv ;
        grid = patch('vertices',[GX(:) GY(:)]...
                        ,'faces', p1(:) + [0 gridDiv gridDiv+1 1] ...
                        ,'facecolor','none','edgecolor','w') ;
        pts = gobjects(0) ;
        for nn = 1:elmt.nNodes
            pts(nn) = drawpoint('Position',elmt.NodeLocalCoordinates(nn,:).*[nJ-1 nI-1] + 1) ;
        end
        callback = @(src,evt)set(grid,'vertices',elmt.evalAt([GX(:) GY(:)])*cat(1,pts.Position)) ;
        for nn = 1:elmt.nNodes
            addlistener(pts(nn),'MovingROI',callback) ;
        end
        callback() ;
        axis equal
        axis tight
        colormap gray
    ax2 = axes('OuterPosition',[.5 0 .5 1]) ;
        srf = surf(JJ,II,JJ*0,refImg,'facecolor','flat','edgecolor','none') ;
        im = imagesc(refImg) ;
        roi = drawrectangle('Position',[1 1 nJ-1 nI-1],'facealpha',0,'linewidth',.5) ;
        %axis ij
        axis equal
        axis tight
        colormap gray
        
%% Transform the second image
% Interpolation matrix
    mesh = pkg.geometry.mesh.Mesh('Nodes',cat(1,pts.Position)...
                ,'Elems',pkg.geometry.mesh.elements.ElementTable(...
                    'Types',elmt,'Indices',[1 1:elmt.nNodes]) ...
                ) ;
    P0 = [JJ(:) II(:)]+.5 ;
    [E,ie] = mesh.localize(P0,mesh.Elems,true) ;
    M = mesh.interpMat(E,ie) ;
% New image
    pp = elmt.evalAt([0 0 ; 1 0 ; 1 1 ; 0 1])*cat(1,pts.Position) ;
    L = sqrt(sum(diff(pp([1:end,1],:),1,1).^2,2)) ;
    L = max(reshape(L,[2 2]),[],2)' ;
    P = M*(elmt.NodeLocalCoordinates.*L) ; 
    P = reshape(P,[nI+1 nJ+1 2]) ;
    srf.XData = P(:,:,1) ;
    srf.YData = P(:,:,2) ;
        
%% Transform the second image
% Create a mesh
    mesh = pkg.geometry.mesh.Mesh('Nodes',cat(1,pts.Position)...
                ,'Elems',pkg.geometry.mesh.elements.ElementTable(...
                    'Types',elmt,'Indices',[1 1:elmt.nNodes]) ...
                ) ;
% Localize the image corners in the mesh
    P0 = [1 1 ; nJ 1 ; nJ nI ; 1 nI] ;
    E = mesh.localize(P0,mesh.Elems,true) ;
    bbox = [min(E,[],1) ; max(E,[],1)] ;
% Create the new grid of points
    % Grid extent
        pp = elmt.evalAt([0 0 ; 1 0 ; 1 1 ; 0 1])*cat(1,pts.Position) ;
        L = sqrt(sum(diff(pp([1:end,1],:),1,1).^2,2)) ;
        L = max(reshape(L,[2 2]),[],2)' ;
    % Grid points
        dx = 1./L ;
        [E1,E2] = meshgrid(bbox(1,1):dx(1):bbox(2,1),bbox(1,2):dx(2):bbox(2,2)) ;
        ji = elmt.evalAt([E1(:) E2(:)])*cat(1,pts.Position) ;
% Interpolate
    newImg = interp2(refImg,ji(:,1),ji(:,2),'linear',0) ;
    newImg = reshape(newImg,size(E1)) ;
% Display
    im.CData = newImg ;
    
%% INTERPOLATE ALL IMAGES
% Crop with roi rectangle
    roiPos = roi.Position(1:2)+[0;1]*roi.Position(3:4) ;
    validE1 = ceil(roiPos(1,1)):floor(roiPos(2,1)) ;
    validE2 = ceil(roiPos(1,2)):floor(roiPos(2,2)) ;
    vE1 = E1(validE2,validE1) ;
    vE2 = E2(validE2,validE1) ;
% Interpolation pixels
    ji = elmt.evalAt([vE1(:) vE2(:)])*cat(1,pts.Position) ;
% New Images
    nIMG = numel(hd.Images{cam}) ;
    IMG = cell(1,nIMG) ;
    wtbr = waitbar(0,'Interpolation...') ;
    for ii = 1:nIMG
        IMG{ii} = interp2(hd.Images{cam}{ii},ji(:,1),ji(:,2),'linear',0) ;
        IMG{ii} = reshape(IMG{ii},size(vE1)) ;
        wtbr = waitbar(ii/nIMG,wtbr) ;
    end
    delete(wtbr) ;
    
%% PUSH TO NAVDIC
    newCamIdx = 3 ; numel(hd.Cameras)+1 ;
    newCam = hd.Cameras(cam) ;
    newCam.Name = ['Transform | ' hd.Cameras(cam).Name] ;
    hd.Images{newCamIdx} = IMG ;
    hd.Cameras(newCamIdx) = newCam ;
    
    
    
    
    
    
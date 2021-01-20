%% CREATE A GRID MESH SEED
    global hd
    cam = 1 ;
    refFrame = hd.nFrames ;
    
    refImg = hd.Images{cam}{refFrame} ;
    [nI,nJ,nC] = size(refImg) ;
        
    clf ;
    imagesc(refImg)
    axis ij
    axis equal
    axis tight
    colormap gray
    
%% Create the grid
    % Number of nodes
        nX = round(1500/25) ; nY = NaN ; round(nI/nJ*(nX-1))+1 ;
    % Initial grid corners position
        pos0 = round( [1 1 ; nJ 1 ; nJ nI ; 1 nI] + 0.02*[1 1 ; -1 1 ; -1 -1 ; 1 -1].*[nJ nI] ) ;
        poly = findobj(gca,'type','images.roi.polygon') ;
        if isempty(poly) 
            poly = images.roi.Polygon('Parent',gca ...
                                    ,'Position',pos0 ...
                                    ,'FaceSelectable',false ...
                                    ,'FaceAlpha',0 ...
                                    ) ; 
        end
    % Grid Function
        if isnan(nY)
            L = range(poly.Position,1) ;
            nY = round(L(2)/L(1)*(nX-1))+1 ;
        end
        [ee,nn] = meshgrid(linspace(0,1,nX),linspace(0,1,nY)) ;
        N = [(1-ee(:)).*(1-nn(:)) ee(:).*(1-nn(:)) ee(:).*nn(:) (1-ee(:)).*nn(:)] ; 
        grid = @(pos) N*pos ;
    % Point plot
        IND = reshape(1:nX*nY,nY,nX) ;
        p1 = IND(1:end-1,1:end-1) ; p2 = p1+1 ; p3 = p2+nY ; p4 = p3-1 ;
        elems = [p1(:) p2(:) p3(:) p4(:)] ;
        pts = findobj(gca,'type','patch') ;
        if isempty(pts) 
            pts = patch(gca,'faces',elems ...
                            ,'vertices',grid(poly.Position) ...
                            ,'hittest','off' ...
                            ,'pickableparts','none' ...
                            ,'facecolor','none' ...
                            ,'edgecolor','b' ...
                            ,'marker','.' ...
                            ,'markeredgecolor','b' ...
                            ,'markersize',25 ...
                            ) ; 
        end
        pts.Faces = elems ;
        pts.Vertices = grid(poly.Position) ;
    % Polygon listener
        list = event.listener(poly,'MovingROI',@(src,evt)set(pts,'vertices',grid(poly.Position))) ;
        
%% Adjust points individually
    delete(findobj(gca,'type','images.roi.point'))
    % Create the control points
        nodes = images.roi.Point.empty ;
        for ii = 1:nX*nY ; nodes(end+1) = images.roi.Point(gca,'Position',pts.Vertices(ii,:)) ; end
    % Listen their position
        for ii = 1:nX*nY ; addlistener(nodes(ii),'MovingROI',@(src,evt)set(pts,'vertices',reshape([nodes.Position],[],nX*nY)')) ; end
        
%% Create the seed
    newSeed = copy(hd.Seeds(end)) ;
    newSeed.Name = 'Grid' ;
    newSeed.Elems = elems ;
    newSeed.Points = pts.Vertices ;
    newSeed.MovingPoints = repmat(pts.Vertices,[1 1 hd.nFrames])*NaN ;
    
%% Push in navDIC
    hd.Seeds(end) = newSeed ;
    
    

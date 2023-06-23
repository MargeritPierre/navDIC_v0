%% CREATE A GRID MESH SEED
    global hd
    cam = 1 ;
    refFrame = hd.CurrentFrame ;
    
    refImg = hd.Images{cam}{refFrame} ;
    [nI,nJ,nC] = size(refImg) ;
        
    clf ;
    imagesc(refImg)
    axis ij
    axis equal
    axis tight
    colormap gray
    
%% REGULAR GRID
    % Number of nodes
        nX = 8 ; round(150/25) ; 
        nY = 15 ; NaN ; round(nI/nJ*(nX-1))+1 ;
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
        
%% HEXAGONAL GRID
    % Number of cells
        nX = 8 ; round(150/25) ; 
        nY = 14 ;
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
        [ee,nn] = meshgrid((0:nX)/nX,(0:nY)/nY) ; % initial corners
        ee = ee + reshape([0 .5 1.5 2]/3/nX,1,1,[]) ;
        nn = nn + reshape([0 .5 .5 0]/nY,1,1,[]) ;
    % Bars
        idx = reshape(1:numel(ee),size(ee)) ;
        bars = [idx(:,:,[1 2]) ; idx(:,:,[2,3]) ; idx(:,:,[3,4])] ;
        bars = reshape(bars,[],2) ;
        bars = [bars ; reshape([idx(:,1:end-1,4) idx(:,2:end,1)],[],2)] ;
        bars = [bars ; reshape([idx(1:end-1,:,2) idx(2:end,:,1)],[],2)] ;
        bars = [bars ; reshape([idx(1:end-1,:,3) idx(2:end,:,4)],[],2)] ;
    % Cull out of frame
        valid = ee(:)<=1+eps & nn(:)<=1+eps ; ee = ee(valid) ; nn = nn(valid) ;
        bars(any(~valid(bars),2),:) = [] ;
        cv = cumsum(valid) ;
        bars = cv(bars) ;
    % Grid warping function
        N = [(1-ee(:)).*(1-nn(:)) ee(:).*(1-nn(:)) ee(:).*nn(:) (1-ee(:)).*nn(:)] ; 
        grid = @(pos) N*pos ;
    % Point plot
        elems = reshape(bars,[],2);%(1:numel(nn))' ;%[p1(:) p2(:) p3(:) p4(:)] ;
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
    newSeed.Name = 'HEX' ;
    newSeed.Elems = elems ;
    newSeed.Points = pts.Vertices ;
    newSeed.MovingPoints = repmat(pts.Vertices,[1 1 hd.nFrames])*NaN ;
    newSeed.CamIDs = cam ;
    newSeed.refImgs = {refImg} ;
    
%% Push in navDIC
    hd.Seeds(end+1) = newSeed ;
    
%% Assemble multiple grids
    seeds = hd.Seeds(3:end) ;
    
    newSeed = copy(seeds(end)) ;
    newSeed.Name = 'AllPoints' ;
    newSeed.Elems = cat(1,seeds.Elems) ;
    newSeed.Points = cat(1,seeds.Points) ;
    newSeed.MovingPoints = cat(1,seeds.MovingPoints) ;
    
    % Delaunay triangulation
    newSeed.Elems = delaunay(newSeed.Points) ;
    xe = reshape(newSeed.Points(newSeed.Elems,:),[],3,2) ;
    newSeed.Elems(polyarea(xe(:,:,1),xe(:,:,2),2)<1,:) = [] ;
    
